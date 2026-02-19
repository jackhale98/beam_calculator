#!/usr/bin/env python3
"""
Continuous Beam Analysis Solver
================================
Solves statically indeterminate continuous beams using the
direct stiffness (finite element) method with Euler-Bernoulli beam theory.

Supports:
  - Multiple support types: pinned, roller, fixed
  - Load types: point loads, distributed loads (uniform/trapezoidal), moments
  - Metric or Imperial unit systems
  - Output: reactions, shear, moment, slope, deflection, shear stress, bending stress

Usage:
  Configure the beam in the CONFIGURATION section at the bottom of this file,
  then run: python beam_solver.py
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


# ─────────────────────────────────────────────────────────────
# Data Classes for Input
# ─────────────────────────────────────────────────────────────

@dataclass
class Support:
    """
    Beam support definition.
    position : distance from left end of beam
    kind     : 'pinned', 'roller', or 'fixed'
               (pinned and roller both restrain vertical displacement;
                fixed also restrains rotation)
    """
    position: float
    kind: str = "pinned"  # 'pinned', 'roller', 'fixed'


@dataclass
class PointLoad:
    """
    Concentrated force.
    position  : distance from left end
    magnitude : force value (positive = DOWNWARD)
    """
    position: float
    magnitude: float


@dataclass
class DistributedLoad:
    """
    Distributed (line) load over a region.
    start   : start position from left end
    end     : end position from left end
    w_start : load intensity at start (positive = DOWNWARD, force/length)
    w_end   : load intensity at end   (positive = DOWNWARD, force/length)
    For uniform load, set w_start = w_end.
    """
    start: float
    end: float
    w_start: float
    w_end: float


@dataclass
class AppliedMoment:
    """
    Concentrated moment.
    position  : distance from left end
    magnitude : moment value (positive = COUNTERCLOCKWISE)
    """
    position: float
    magnitude: float


@dataclass
class BeamSection:
    """
    Cross-section and material properties.
    E : modulus of elasticity
    I : second moment of area (moment of inertia)
    A : cross-sectional area
    c : distance from neutral axis to extreme fiber (for bending stress)
    """
    E: float
    I: float
    A: float
    c: float


# ─────────────────────────────────────────────────────────────
# Solver Class
# ─────────────────────────────────────────────────────────────

class ContinuousBeamSolver:
    """
    Finite element solver for continuous beams.
    Uses Euler-Bernoulli beam elements with 2 DOFs per node
    (vertical displacement and rotation).
    """

    def __init__(
        self,
        length: float,
        section: BeamSection,
        supports: List[Support],
        point_loads: List[PointLoad] = None,
        distributed_loads: List[DistributedLoad] = None,
        applied_moments: List[AppliedMoment] = None,
        unit_system: str = "imperial",
        n_elements: int = 800,
    ):
        self.length = length
        self.section = section
        self.supports = supports
        self.point_loads = point_loads or []
        self.distributed_loads = distributed_loads or []
        self.applied_moments = applied_moments or []
        self.unit_system = unit_system.lower()
        self.n_elements = n_elements

        # Results (populated after solve)
        self.nodes = None
        self.U = None
        self.reactions = None
        self.x_plot = None
        self.V_plot = None
        self.M_plot = None
        self.theta_plot = None
        self.delta_plot = None
        self.tau_plot = None
        self.sigma_plot = None

    # ── Mesh Generation ──────────────────────────────────────

    def _build_mesh(self) -> np.ndarray:
        """Create a fine mesh with nodes at all key points."""
        key_points = set([0.0, self.length])
        for s in self.supports:
            key_points.add(s.position)
        for p in self.point_loads:
            key_points.add(p.position)
        for d in self.distributed_loads:
            key_points.add(d.start)
            key_points.add(d.end)
        for m in self.applied_moments:
            key_points.add(m.position)

        key_points = sorted(key_points)
        all_nodes = []

        for i in range(len(key_points) - 1):
            x0, x1 = key_points[i], key_points[i + 1]
            seg_len = x1 - x0
            if seg_len < 1e-14:
                continue
            n_seg = max(int(self.n_elements * seg_len / self.length), 4)
            seg_nodes = np.linspace(x0, x1, n_seg + 1)
            if len(all_nodes) > 0:
                seg_nodes = seg_nodes[1:]  # avoid duplicate
            all_nodes.extend(seg_nodes.tolist())

        return np.array(all_nodes)

    # ── Assembly & Solution ──────────────────────────────────

    def solve(self):
        """Assemble and solve the FEM system, then post-process."""
        nodes = self._build_mesh()
        n_nodes = len(nodes)
        n_dof = 2 * n_nodes
        EI = self.section.E * self.section.I

        # ── Global stiffness matrix and force vector ──
        K = np.zeros((n_dof, n_dof))
        F = np.zeros(n_dof)

        for i in range(n_nodes - 1):
            x0 = nodes[i]
            x1 = nodes[i + 1]
            h = x1 - x0
            if h < 1e-15:
                continue

            # Element stiffness matrix
            k_e = (EI / h ** 3) * np.array(
                [
                    [12, 6 * h, -12, 6 * h],
                    [6 * h, 4 * h ** 2, -6 * h, 2 * h ** 2],
                    [-12, -6 * h, 12, -6 * h],
                    [6 * h, 2 * h ** 2, -6 * h, 4 * h ** 2],
                ]
            )

            dofs = [2 * i, 2 * i + 1, 2 * (i + 1), 2 * (i + 1) + 1]
            for a in range(4):
                for b in range(4):
                    K[dofs[a], dofs[b]] += k_e[a, b]

            # ── Distributed loads (equivalent nodal forces) ──
            for dl in self.distributed_loads:
                if x1 <= dl.start + 1e-12 or x0 >= dl.end - 1e-12:
                    continue

                # Clamp element to load region
                xa = max(x0, dl.start)
                xb = min(x1, dl.end)
                load_span = dl.end - dl.start if abs(dl.end - dl.start) > 1e-14 else 1.0

                # Linear interpolation of load at element boundaries
                frac_a = (xa - dl.start) / load_span
                frac_b = (xb - dl.start) / load_span
                wa = dl.w_start + frac_a * (dl.w_end - dl.w_start)
                wb = dl.w_start + frac_b * (dl.w_end - dl.w_start)
                w_avg = (wa + wb) / 2.0
                w_len = xb - xa

                # Approximate as UDL over this sub-element
                # (accurate for fine mesh)
                # Convention: positive load = downward = negative force
                # in upward-positive DOF convention
                # Place load as if acting over the full element length h
                # but only on the portion from xa to xb
                # For simplicity with fine mesh, use lumped loads
                f_left = -w_avg * w_len / 2.0
                f_right = -w_avg * w_len / 2.0
                m_left = -w_avg * w_len ** 2 / 12.0
                m_right = w_avg * w_len ** 2 / 12.0

                # Adjust if load doesn't cover full element
                if abs(w_len - h) > 1e-12:
                    # Use simple lumped approach for partial coverage
                    center = (xa + xb) / 2.0
                    rel = (center - x0) / h  # 0 to 1
                    total_force = -w_avg * w_len
                    f_left = total_force * (1.0 - rel)
                    f_right = total_force * rel
                    m_left = 0.0
                    m_right = 0.0

                F[dofs[0]] += f_left
                F[dofs[1]] += m_left
                F[dofs[2]] += f_right
                F[dofs[3]] += m_right

        # ── Point loads ──
        for pl in self.point_loads:
            idx = np.argmin(np.abs(nodes - pl.position))
            F[2 * idx] -= pl.magnitude  # positive load = downward = negative DOF force

        # ── Applied moments ──
        for m in self.applied_moments:
            idx = np.argmin(np.abs(nodes - m.position))
            F[2 * idx + 1] += m.magnitude  # positive CCW moment

        # ── Boundary conditions ──
        constrained_dofs = []
        support_node_indices = {}
        for s in self.supports:
            idx = np.argmin(np.abs(nodes - s.position))
            support_node_indices[s.position] = idx
            constrained_dofs.append(2 * idx)  # v = 0
            if s.kind == "fixed":
                constrained_dofs.append(2 * idx + 1)  # theta = 0

        constrained_dofs = sorted(set(constrained_dofs))
        free_dofs = [i for i in range(n_dof) if i not in constrained_dofs]

        # ── Solve ──
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        F_f = F[free_dofs]

        try:
            U = np.zeros(n_dof)
            U[free_dofs] = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                "Stiffness matrix is singular. Check that the beam is "
                "properly supported (not a mechanism)."
            )

        # ── Reactions ──
        R = K @ U - F
        reactions = {}
        for s in self.supports:
            idx = support_node_indices[s.position]
            vert = R[2 * idx]
            mom = R[2 * idx + 1] if s.kind == "fixed" else 0.0
            reactions[s.position] = {"vertical": vert, "moment": mom}

        # ── Store raw results ──
        self.nodes = nodes
        self.U = U
        self.reactions = reactions

        # ── Post-process: V, M, theta, delta, stress ──
        self._post_process()

        return self

    # ── Post-Processing ──────────────────────────────────────

    def _distributed_load_at(self, x: float) -> float:
        """Return the total distributed load intensity at position x (positive downward)."""
        w_total = 0.0
        for dl in self.distributed_loads:
            if dl.start - 1e-12 <= x <= dl.end + 1e-12:
                span = dl.end - dl.start
                if span < 1e-14:
                    continue
                frac = (x - dl.start) / span
                frac = max(0.0, min(1.0, frac))
                w_total += dl.w_start + frac * (dl.w_end - dl.w_start)
        return w_total

    def _post_process(self):
        """Compute shear, moment, slope, deflection, and stress arrays.

        Uses Hermite shape function derivatives to extract M and V
        directly from the FEM displacement field, which is more robust
        than an equilibrium walk.

        For Euler-Bernoulli beam elements with DOFs [v1, θ1, v2, θ2]:
          v(ξ) = N1·v1 + N2·θ1 + N3·v2 + N4·θ2   (ξ = (x-x1)/h)

          M = EI·v'' :
            v''(ξ) = (1/h²)[(-6+12ξ)v1 + h(-4+6ξ)θ1 + (6-12ξ)v2 + h(-2+6ξ)θ2]

          V = dM/dx = EI·v''' :
            v'''(ξ) = (1/h³)[12·v1 + 6h·θ1 - 12·v2 + 6h·θ2]
            (constant within each element)
        """
        nodes = self.nodes
        U = self.U
        n_pts = len(nodes)
        EI = self.section.E * self.section.I

        # ── Deflection and slope from FEM solution ──
        delta = np.array([U[2 * i] for i in range(n_pts)])
        theta = np.array([U[2 * i + 1] for i in range(n_pts)])

        V = np.zeros(n_pts)
        M = np.zeros(n_pts)

        # ── Compute M and V from shape function derivatives ──
        # At each node, evaluate using the element to its RIGHT.
        # This gives values "just to the right" of each node,
        # naturally capturing discontinuities at supports and loads.
        for i in range(n_pts - 1):
            h = nodes[i + 1] - nodes[i]
            if h < 1e-15:
                continue
            v1 = U[2 * i]
            t1 = U[2 * i + 1]
            v2 = U[2 * (i + 1)]
            t2 = U[2 * (i + 1) + 1]

            # M at ξ=0 (left end of element)
            M[i] = (EI / h ** 2) * (-6 * v1 - 4 * h * t1 + 6 * v2 - 2 * h * t2)

            # V (constant within element)
            V[i] = (EI / h ** 3) * (12 * v1 + 6 * h * t1 - 12 * v2 + 6 * h * t2)

        # Last node: evaluate using element to its LEFT at ξ=1
        if n_pts >= 2:
            i = n_pts - 2
            h = nodes[i + 1] - nodes[i]
            if h > 1e-15:
                v1 = U[2 * i]
                t1 = U[2 * i + 1]
                v2 = U[2 * (i + 1)]
                t2 = U[2 * (i + 1) + 1]
                M[-1] = (EI / h ** 2) * (6 * v1 + 2 * h * t1 - 6 * v2 + 4 * h * t2)
                V[-1] = (EI / h ** 3) * (12 * v1 + 6 * h * t1 - 12 * v2 + 6 * h * t2)

        # ── Stresses ──
        tau = V / self.section.A  # average shear stress
        sigma = M * self.section.c / self.section.I  # max bending stress (at extreme fiber)

        # ── Store ──
        self.x_plot = nodes
        self.V_plot = V
        self.M_plot = M
        self.theta_plot = theta
        self.delta_plot = delta
        self.tau_plot = tau
        self.sigma_plot = sigma

    # ── Printing ──────────────────────────────────────────────

    def print_results(self):
        """Print reactions and key values to console."""
        if self.unit_system == "imperial":
            f_unit = "lb"
            l_unit = "in"
            m_unit = "lb·in"
            s_unit = "psi"
            d_unit = "in"
            ang_unit = "rad"
        else:
            f_unit = "N"
            l_unit = "mm"
            m_unit = "N·mm"
            s_unit = "MPa"
            d_unit = "mm"
            ang_unit = "rad"

        print("=" * 65)
        print("  CONTINUOUS BEAM ANALYSIS RESULTS")
        print(f"  Unit system: {self.unit_system}")
        print("=" * 65)

        print("\n── SUPPORT REACTIONS ─────────────────────────────────────")
        total_v = 0.0
        for pos in sorted(self.reactions.keys()):
            r = self.reactions[pos]
            v = r["vertical"]
            m = r["moment"]
            total_v += v
            support_type = "unknown"
            for s in self.supports:
                if abs(s.position - pos) < 1e-10:
                    support_type = s.kind
                    break
            print(f"  x = {pos:10.4f} {l_unit}  ({support_type})")
            print(f"      Vertical reaction: {v:12.4f} {f_unit}  ({'↑' if v > 0 else '↓'})")
            if abs(m) > 1e-10:
                print(f"      Moment reaction:   {m:12.4f} {m_unit}  ({'↺' if m > 0 else '↻'})")

        print(f"\n  Sum of vertical reactions: {total_v:.4f} {f_unit}")

        # Total applied load
        total_load = sum(pl.magnitude for pl in self.point_loads)
        for dl in self.distributed_loads:
            w_avg = (dl.w_start + dl.w_end) / 2
            total_load += w_avg * (dl.end - dl.start)
        print(f"  Total applied load:       {total_load:.4f} {f_unit} (downward)")
        print(f"  Equilibrium check:        {abs(total_v - total_load):.6e} {f_unit}")

        print("\n── EXTREME VALUES ───────────────────────────────────────")
        print(f"  Max shear force:     {np.max(np.abs(self.V_plot)):12.4f} {f_unit}")
        print(f"  Max bending moment:  {np.max(np.abs(self.M_plot)):12.4f} {m_unit}")
        print(f"  Max deflection:      {np.max(np.abs(self.delta_plot)):12.6f} {d_unit}")
        print(f"  Max slope:           {np.max(np.abs(self.theta_plot)):12.6e} {ang_unit}")
        print(f"  Max avg shear stress:{np.max(np.abs(self.tau_plot)):12.4f} {s_unit}")
        print(f"  Max bending stress:  {np.max(np.abs(self.sigma_plot)):12.4f} {s_unit}")
        print("=" * 65)

    # ── Plotting ──────────────────────────────────────────────

    def plot(self, save_path: Optional[str] = None):
        """Generate all diagrams."""
        if self.unit_system == "imperial":
            f_unit = "lb"
            l_unit = "in"
            m_unit = "lb·in"
            s_unit = "psi"
            d_unit = "in"
        else:
            f_unit = "N"
            l_unit = "mm"
            m_unit = "N·mm"
            s_unit = "MPa"
            d_unit = "mm"

        x = self.x_plot

        fig, axes = plt.subplots(6, 1, figsize=(14, 24), sharex=True)
        fig.suptitle("Continuous Beam Analysis", fontsize=16, fontweight="bold", y=0.995)

        # ── Helper: draw beam schematic on each subplot ──
        def draw_supports(ax):
            y_min, y_max = ax.get_ylim()
            for s in self.supports:
                ax.axvline(s.position, color="gray", linewidth=0.5, linestyle="--", alpha=0.4)

        plot_data = [
            (self.V_plot, f"Shear Force ({f_unit})", "steelblue", "Shear Force Diagram"),
            (self.M_plot, f"Bending Moment ({m_unit})", "firebrick", "Bending Moment Diagram"),
            (self.theta_plot, f"Slope (rad)", "darkorange", "Slope Diagram"),
            (self.delta_plot, f"Deflection ({d_unit})", "seagreen", "Deflection Diagram"),
            (self.tau_plot, f"Avg Shear Stress ({s_unit})", "mediumpurple", "Average Shear Stress Diagram"),
            (self.sigma_plot, f"Bending Stress ({s_unit})", "crimson", "Bending Stress Diagram"),
        ]

        for idx, (data, ylabel, color, title) in enumerate(plot_data):
            ax = axes[idx]
            ax.fill_between(x, 0, data, alpha=0.25, color=color)
            ax.plot(x, data, color=color, linewidth=1.5)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_ylabel(ylabel, fontsize=10)
            ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
            ax.grid(True, alpha=0.3)

            # Draw support markers
            for s in self.supports:
                ax.axvline(s.position, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)

            # Mark extremes
            i_max = np.argmax(np.abs(data))
            ax.annotate(
                f"{data[i_max]:.4g}",
                xy=(x[i_max], data[i_max]),
                fontsize=8,
                color=color,
                fontweight="bold",
                ha="center",
                va="bottom" if data[i_max] >= 0 else "top",
            )

        axes[-1].set_xlabel(f"Position along beam ({l_unit})", fontsize=11)

        # ── Draw beam schematic at the top ──
        ax_beam = fig.add_axes([0.125, 0.86, 0.775, 0.12])
        beam_len = self.length
        pad = beam_len * 0.03
        ax_beam.set_xlim(x[0] - pad, x[-1] + pad)
        # Coordinate layout:
        #   beam line at y = 0
        #   loads above (positive y)
        #   supports below (negative y), reactions below supports
        ax_beam.set_ylim(-6.5, 5.5)
        ax_beam.axhline(0, color="black", linewidth=3, zorder=3)
        ax_beam.set_yticks([])
        ax_beam.set_xticks([])
        ax_beam.set_frame_on(False)
        ax_beam.set_title("Beam Schematic", fontsize=10, fontweight="bold", loc="left")

        # ── Supports (drawn BELOW the beam line) with reaction labels ──
        for s in self.supports:
            sx = s.position
            r = self.reactions.get(sx, {"vertical": 0, "moment": 0})
            rv = r["vertical"]
            rm = r["moment"]

            if s.kind == "fixed":
                # Fixed support: filled square below beam + hatching
                rect_w = beam_len * 0.02
                ax_beam.add_patch(plt.Rectangle(
                    (sx - rect_w, -1.8), 2 * rect_w, 1.8,
                    facecolor="dimgray", edgecolor="black", linewidth=1.2, zorder=4
                ))
                # Hatching lines
                for yh in np.linspace(-1.8, -0.3, 5):
                    ax_beam.plot(
                        [sx - rect_w * 1.5, sx + rect_w * 1.5],
                        [yh - 0.3, yh],
                        color="black", linewidth=0.7, zorder=4
                    )
                support_bottom = -1.8
                label_base_y = -2.5
            elif s.kind == "roller":
                # Roller: triangle + circle below
                tri_h = 1.2
                tri_w = beam_len * 0.025
                triangle = plt.Polygon(
                    [[sx, 0], [sx - tri_w, -tri_h], [sx + tri_w, -tri_h]],
                    facecolor="white", edgecolor="black", linewidth=1.2, zorder=4
                )
                ax_beam.add_patch(triangle)
                circle = plt.Circle((sx, -tri_h - 0.3), 0.25,
                                    facecolor="white", edgecolor="black", linewidth=1.0, zorder=4)
                ax_beam.add_patch(circle)
                ax_beam.plot([sx - tri_w, sx + tri_w], [-tri_h - 0.6, -tri_h - 0.6],
                             color="black", linewidth=1.0, zorder=4)
                support_bottom = -tri_h - 0.6
                label_base_y = -2.3
            else:
                # Pinned: triangle only
                tri_h = 1.2
                tri_w = beam_len * 0.025
                triangle = plt.Polygon(
                    [[sx, 0], [sx - tri_w, -tri_h], [sx + tri_w, -tri_h]],
                    facecolor="white", edgecolor="black", linewidth=1.2, zorder=4
                )
                ax_beam.add_patch(triangle)
                ax_beam.plot([sx - tri_w * 1.2, sx + tri_w * 1.2], [-tri_h, -tri_h],
                             color="black", linewidth=1.2, zorder=4)
                support_bottom = -tri_h
                label_base_y = -2.3

            # ── Reaction force arrow and label ──
            arrow_top = support_bottom - 0.2
            arrow_bot = arrow_top - 1.8
            direction = "↑" if rv > 0 else "↓"
            arrow_color = "darkgreen" if rv > 0 else "orangered"

            ax_beam.annotate(
                "",
                xy=(sx, arrow_top),
                xytext=(sx, arrow_bot),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.15",
                                color=arrow_color, lw=2),
                zorder=5,
            )

            # Position and type label
            label_parts = [f"{sx:.4g} {l_unit} ({s.kind})"]
            # Reaction value label
            label_parts.append(f"R = {rv:.2f} {f_unit} {direction}")
            if s.kind == "fixed" and abs(rm) > 1e-10:
                m_dir = "↺" if rm > 0 else "↻"
                label_parts.append(f"M = {rm:.2f} {m_unit} {m_dir}")

            ax_beam.text(
                sx, arrow_bot - 0.3,
                "\n".join(label_parts),
                ha="center", va="top", fontsize=6.5, color=arrow_color,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=arrow_color, alpha=0.9),
                zorder=6,
            )

        # ── Distributed loads (drawn ABOVE the beam line) ──
        for dl in self.distributed_loads:
            load_span = dl.end - dl.start
            if load_span < 1e-14:
                continue
            w_max_abs = max(abs(dl.w_start), abs(dl.w_end), 1e-14)
            arrow_scale = 2.5  # max arrow height in axes coords

            # Draw arrows
            n_arrows = max(int(10 * load_span / beam_len), 5)
            x_arr = np.linspace(dl.start, dl.end, n_arrows)
            for xa in x_arr:
                frac = (xa - dl.start) / load_span
                w_here = dl.w_start + frac * (dl.w_end - dl.w_start)
                arr_h = arrow_scale * abs(w_here) / w_max_abs
                arr_h = max(arr_h, 0.3)
                ax_beam.annotate(
                    "",
                    xy=(xa, 0.15),
                    xytext=(xa, 0.15 + arr_h),
                    arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.15",
                                    color="royalblue", lw=1.2),
                    zorder=5,
                )

            # Draw connecting line across arrow tops
            x_top = np.linspace(dl.start, dl.end, 50)
            y_top = []
            for xt in x_top:
                frac = (xt - dl.start) / load_span
                w_h = dl.w_start + frac * (dl.w_end - dl.w_start)
                y_top.append(0.15 + arrow_scale * abs(w_h) / w_max_abs)
            ax_beam.plot(x_top, y_top, color="royalblue", linewidth=1.5, zorder=5)

            # Label
            if abs(dl.w_start - dl.w_end) < 1e-10:
                label = f"w = {dl.w_start:.4g} {f_unit}/{l_unit}\n@ {dl.start:.4g}–{dl.end:.4g} {l_unit}"
            else:
                label = f"w = {dl.w_start:.4g}→{dl.w_end:.4g} {f_unit}/{l_unit}\n@ {dl.start:.4g}–{dl.end:.4g} {l_unit}"
            mid_x = (dl.start + dl.end) / 2
            frac_mid = 0.5
            w_mid = dl.w_start + frac_mid * (dl.w_end - dl.w_start)
            label_y = 0.15 + arrow_scale * abs(w_mid) / w_max_abs + 0.5
            ax_beam.text(mid_x, label_y, label,
                         ha="center", va="bottom", fontsize=8, color="royalblue",
                         fontweight="bold",
                         bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                   edgecolor="royalblue", alpha=0.9))

        # ── Point loads (drawn ABOVE the beam line) ──
        for pl in self.point_loads:
            arr_h = 2.5
            ax_beam.annotate(
                "",
                xy=(pl.position, 0.15),
                xytext=(pl.position, 0.15 + arr_h),
                arrowprops=dict(arrowstyle="->,head_width=0.4,head_length=0.2",
                                color="red", lw=2.5),
                zorder=6,
            )
            direction = "↓" if pl.magnitude > 0 else "↑"
            ax_beam.text(
                pl.position, 0.15 + arr_h + 0.3,
                f"P = {abs(pl.magnitude):.4g} {f_unit} {direction}\n@ {pl.position:.4g} {l_unit}",
                ha="center", va="bottom", fontsize=8, color="red", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="red", alpha=0.9),
                zorder=6,
            )

        # ── Applied moments (drawn ABOVE the beam line) ──
        for am in self.applied_moments:
            # Draw curved arrow for moment
            arc_radius = beam_len * 0.02
            if am.magnitude > 0:
                arc_style = "arc3,rad=-0.5"
                direction = "↺ CCW"
            else:
                arc_style = "arc3,rad=0.5"
                direction = "↻ CW"

            # Draw a circular arc arrow
            theta_vals = np.linspace(0, 1.5 * np.pi, 40)
            cx, cy = am.position, 1.8
            r = beam_len * 0.018
            arc_x = cx + r * np.cos(theta_vals)
            arc_y = cy + r * np.sin(theta_vals)
            ax_beam.plot(arc_x, arc_y, color="green", linewidth=2, zorder=6)
            # Arrowhead at end of arc
            dx = arc_x[-1] - arc_x[-2]
            dy = arc_y[-1] - arc_y[-2]
            ax_beam.annotate(
                "",
                xy=(arc_x[-1], arc_y[-1]),
                xytext=(arc_x[-3], arc_y[-3]),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.15",
                                color="green", lw=2),
                zorder=6,
            )
            ax_beam.text(
                am.position, cy + r + 0.6,
                f"M = {abs(am.magnitude):.4g} {m_unit} {direction}\n@ {am.position:.4g} {l_unit}",
                ha="center", va="bottom", fontsize=8, color="green", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="green", alpha=0.9),
                zorder=6,
            )

        plt.tight_layout(rect=[0, 0, 1, 0.86])

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"\nPlot saved to: {save_path}")
        else:
            plt.show()

        return fig


# ─────────────────────────────────────────────────────────────
# CONFIGURATION — EDIT THIS SECTION FOR YOUR BEAM
# ─────────────────────────────────────────────────────────────

def example_four_support_beam():
    """
    Example: Continuous beam with 4 supports and a centered distributed load.

    This example uses IMPERIAL units.
    Lengths in inches, forces in pounds, stresses in psi.

    Adjust all values below to match your specific problem.
    """

    # ── Unit system ──
    # 'imperial' : lengths in inches, forces in lb, E in psi, stress in psi
    # 'metric'   : lengths in mm, forces in N, E in MPa, stress in MPa
    unit_system = "imperial"

    # ── Beam total length ──
    beam_length = 120.0  # inches (10 feet)

    # ── Cross-section and material ──
    # Example: Steel W8x10 (approximate properties)
    #   E = 29,000,000 psi (steel)
    #   I = 30.8 in^4
    #   A = 2.96 in^2
    #   c = 3.94 in (half-depth)
    section = BeamSection(
        E=29_000_000.0,  # psi
        I=30.8,          # in^4
        A=2.96,          # in^2
        c=3.94,          # in
    )

    # ── Supports ──
    # Types: 'pinned', 'roller', 'fixed'
    # Pinned and roller both prevent vertical displacement.
    # Fixed also prevents rotation.
    # NOTE: You need at least enough supports to prevent rigid-body motion.
    #        Typically: one pinned + rollers, or one fixed end, etc.
    supports = [
        Support(position=0.0,   kind="pinned"),   # Support A
        Support(position=36.0,  kind="roller"),    # Support B (3 ft)
        Support(position=84.0,  kind="roller"),    # Support C (7 ft)
        Support(position=120.0, kind="roller"),    # Support D (10 ft)
    ]

    # ── Loads ──
    # Point loads: positive magnitude = DOWNWARD
    point_loads = [
        PointLoad(position=18.0, magnitude=200.0),  # 200 lb at 18 in
    ]

    # Distributed loads: positive w = DOWNWARD (force per unit length)
    # For uniform load, set w_start = w_end
    # For trapezoidal, set different values
    distributed_loads = [
        DistributedLoad(start=36.0, end=84.0, w_start=12.5, w_end=12.5),
        # This is 12.5 lb/in over 48 inches = 600 lb total
    ]

    # Applied moments: positive = counterclockwise
    applied_moments = [
        AppliedMoment(position=102.0, magnitude=1500.0),  # 1500 lb·in CCW
    ]

    # ── Create solver and run ──
    solver = ContinuousBeamSolver(
        length=beam_length,
        section=section,
        supports=supports,
        point_loads=point_loads,
        distributed_loads=distributed_loads,
        applied_moments=applied_moments,
        unit_system=unit_system,
        n_elements=800,
    )

    solver.solve()
    solver.print_results()
    solver.plot(save_path="beam_analysis.png")

    return solver


def example_metric_beam():
    """
    Example: Metric units (SI).
    Lengths in mm, forces in N, E in MPa, stress in MPa.
    """
    section = BeamSection(
        E=200_000.0,   # MPa (steel)
        I=8_356_000.0, # mm^4  (IPE 200 approx)
        A=2_850.0,     # mm^2
        c=100.0,       # mm (half-depth)
    )

    supports = [
        Support(position=0.0,    kind="fixed"),
        Support(position=2000.0, kind="roller"),
        Support(position=4000.0, kind="roller"),
        Support(position=6000.0, kind="pinned"),
    ]

    distributed_loads = [
        DistributedLoad(start=1000.0, end=5000.0, w_start=10.0, w_end=10.0),
        # 10 N/mm = 10 kN/m over 4 m = 40 kN total
    ]

    solver = ContinuousBeamSolver(
        length=6000.0,
        section=section,
        supports=supports,
        distributed_loads=distributed_loads,
        unit_system="metric",
        n_elements=800,
    )

    solver.solve()
    solver.print_results()
    solver.plot(save_path="beam_analysis_metric.png")

    return solver


# ─────────────────────────────────────────────────────────────
# MAIN — Run the example
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n" + "=" * 65)
    print("  RUNNING IMPERIAL EXAMPLE (4-support beam, centered UDL)")
    print("=" * 65 + "\n")
    solver = example_four_support_beam()

    # Uncomment the line below to also run the metric example:
    # solver_metric = example_metric_beam()
