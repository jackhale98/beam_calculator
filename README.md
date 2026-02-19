# Continuous Beam Analysis Solver

A Python tool for analyzing statically indeterminate continuous beams using the finite element method (Euler-Bernoulli beam theory). It computes support reactions, internal forces, deflections, and stresses, and produces a complete set of annotated engineering diagrams.

![Example Output](beam_analysis.png)

---

## Features

- **Any number of supports** at arbitrary positions — pinned, roller, fixed, or elastic (beam/spring)
- **Mixed loading** — point loads, distributed loads (uniform or trapezoidal), and applied moments, in any combination
- **Imperial or metric** unit systems
- **Automatic self-weight** — optional beam weight per unit length is automatically included as a uniform distributed load
- **Maximum shear stress** — reports both average (V/A) and maximum shear stress using a configurable cross-section shape factor
- **Seven output diagrams** — shear force, bending moment, slope, deflection, average shear stress, maximum shear stress, and bending stress
- **Annotated beam schematic** — a free body diagram showing all applied loads with positions, support types, and calculated reaction forces with direction arrows
- **Console output** — tabulated reactions, equilibrium check, and extreme values

---

## Requirements

- Python 3.8+
- NumPy
- Matplotlib

Install dependencies:

```bash
pip install numpy matplotlib
```

---

## Quick Start

1. Open `beam_solver.py` and edit the configuration in one of the example functions at the bottom of the file.
2. Run:

```bash
python beam_solver.py
```

3. View the console output for reactions and extreme values, and the saved `beam_analysis.png` for diagrams.

By default, the script runs `example_four_support_beam()`. To run a different example, uncomment the corresponding line in the `if __name__ == "__main__"` block:

```python
if __name__ == "__main__":
    solver = example_four_support_beam()       # 4-support beam with UDL (imperial)

    # Uncomment the line below to also run the metric example:
    # solver_metric = example_metric_beam()

    # Uncomment the line below to run the beam (spring) support example:
    # solver_spring = example_beam_support()
```

| Example | Description |
|---------|-------------|
| `example_four_support_beam()` | Four-support continuous beam with a centered uniform distributed load (imperial units) |
| `example_metric_beam()` | Beam analysis using metric units (N, mm, MPa) |
| `example_beam_support()` | Beam with an elastic spring support at mid-span (`kind="beam"`, `stiffness=5000 lb/in`) |

---

## How to Configure Your Beam

All configuration is done by editing the example function near the bottom of the script. The following sections walk through each input.

### Unit System

```python
unit_system = "imperial"   # lengths: in, forces: lb, E: psi, stress: psi
unit_system = "metric"     # lengths: mm, forces: N,  E: MPa, stress: MPa
```

All inputs and outputs must use consistent units within the chosen system.

### Beam Length

```python
beam_length = 120.0  # total span in inches (or mm for metric)
```

### Cross-Section and Material Properties

```python
section = BeamSection(
    E=29_000_000.0,       # modulus of elasticity (psi or MPa)
    I=30.8,               # moment of inertia (in⁴ or mm⁴)
    A=2.96,               # cross-sectional area (in² or mm²)
    c=3.94,               # distance from neutral axis to extreme fiber (in or mm)
    weight_per_length=0.833,  # self-weight per unit length (lb/in or N/mm), default 0
    shear_shape_factor=1.5,   # multiplier for max shear stress (default 1.5)
)
```

| Property | Description | Imperial | Metric | Default |
|----------|-------------|----------|--------|---------|
| `E` | Modulus of elasticity | psi | MPa | — |
| `I` | Second moment of area | in⁴ | mm⁴ | — |
| `A` | Cross-sectional area | in² | mm² | — |
| `c` | Extreme fiber distance | in | mm | — |
| `weight_per_length` | Beam self-weight per unit length | lb/in | N/mm | 0 |
| `shear_shape_factor` | Max-to-average shear stress ratio | — | — | 1.5 |

#### Self-Weight

If `weight_per_length` is set to a value greater than zero, the solver automatically adds a full-span uniform downward distributed load equal to the beam's self-weight. This load participates in the analysis alongside all other applied loads and is included in the equilibrium check.

```python
# W8x10 steel beam weighs 10 lb/ft = 0.833 lb/in
section = BeamSection(E=29e6, I=30.8, A=2.96, c=3.94, weight_per_length=0.833)
```

To omit self-weight (the default), either leave the parameter out or set it to `0`.

#### Shear Shape Factor

The average shear stress `V/A` underestimates the actual maximum shear stress at the neutral axis. The `shear_shape_factor` converts average to maximum: `τ_max = shear_shape_factor × V/A`.

| Cross-Section | Shape Factor |
|---------------|-------------|
| Rectangular | 1.5 |
| Solid circular | 4/3 ≈ 1.333 |
| Thin-walled circular | 2.0 |
| Wide-flange (I-beam) | ≈ A / (d × t_w) |
| Aluminum extrusion (e.g. 80/20) | ≈ A / (d × t_w_min) — use the thinnest web |

For complex profiles like 80/20 T-slot extrusions, check the manufacturer's data sheet for section properties and minimum web thickness, then compute the factor as `A / (depth × min_web_thickness)`.

**Common material values:**

| Material | E (Imperial) | E (Metric) |
|----------|-------------|------------|
| Steel | 29,000,000 psi | 200,000 MPa |
| Aluminum | 10,000,000 psi | 69,000 MPa |
| Douglas Fir | 1,700,000 psi | 12,000 MPa |

### Supports

```python
supports = [
    Support(position=0.0,   kind="pinned"),
    Support(position=36.0,  kind="roller"),
    Support(position=84.0,  kind="roller"),
    Support(position=120.0, kind="roller"),
]
```

| Type | Vertical Displacement | Rotation | Typical Use |
|------|----------------------|----------|-------------|
| `"pinned"` | Restrained | Free | One per beam (prevents horizontal drift) |
| `"roller"` | Restrained | Free | Additional interior/end supports |
| `"fixed"` | Restrained | Restrained | Cantilever root, wall embedment |
| `"beam"` | Elastic (spring) | Free | Flexible supports, elastic foundations |

#### Beam (Spring) Supports

A `"beam"` support models an elastic/spring support with finite stiffness. Unlike pinned or roller supports which enforce zero displacement, a beam support allows the beam to deflect at that location. The reaction force is proportional to the displacement: `R = k * delta`.

The `stiffness` parameter (in force per unit displacement, e.g. lb/in or N/mm) is required for beam supports and must be positive.

```python
supports = [
    Support(position=0.0,   kind="pinned"),
    Support(position=60.0,  kind="beam", stiffness=5000.0),  # 5000 lb/in spring
    Support(position=120.0, kind="roller"),
]
```

The solver adds the spring stiffness to the global stiffness matrix and solves for the displacement simultaneously with all other DOFs. This means the deflection at the beam support accounts for the beam's own bending stiffness, the spring stiffness, all applied loads, and all other supports — not just a simple `F/k` calculation.

**Stability requirement:** The beam must have enough supports to prevent rigid-body motion. Typically this means at least one pinned support plus additional rollers, or at least one fixed support. Beam (spring) supports alone do not prevent rigid-body motion — they must be combined with at least one pinned or fixed support.

### Point Loads

Concentrated forces applied at a single location. Positive magnitude = **downward**.

```python
point_loads = [
    PointLoad(position=18.0, magnitude=200.0),   # 200 lb downward at 18 in
    PointLoad(position=60.0, magnitude=-100.0),   # 100 lb upward at 60 in
]
```

### Distributed Loads

Line loads applied over a region. Positive intensity = **downward**, in force per unit length.

```python
# Uniform load: set w_start = w_end
distributed_loads = [
    DistributedLoad(start=36.0, end=84.0, w_start=12.5, w_end=12.5),
]

# Trapezoidal load: set different start/end intensities
distributed_loads = [
    DistributedLoad(start=0.0, end=60.0, w_start=5.0, w_end=20.0),
]

# Triangular load: set one end to zero
distributed_loads = [
    DistributedLoad(start=0.0, end=48.0, w_start=0.0, w_end=15.0),
]
```

Multiple distributed loads can overlap — their effects are superimposed.

### Applied Moments

Concentrated moments applied at a single location. Positive magnitude = **counterclockwise**.

```python
applied_moments = [
    AppliedMoment(position=102.0, magnitude=1500.0),   # 1500 lb·in CCW
    AppliedMoment(position=60.0,  magnitude=-800.0),    # 800 lb·in CW
]
```

### Solver Settings

```python
solver = ContinuousBeamSolver(
    length=beam_length,
    section=section,
    supports=supports,
    point_loads=point_loads,
    distributed_loads=distributed_loads,
    applied_moments=applied_moments,
    unit_system=unit_system,
    n_elements=800,  # mesh density — increase for smoother plots
)
```

The `n_elements` parameter controls mesh density. Higher values give smoother plots and marginally better accuracy at discontinuities, at the cost of computation time. 800 is a good default; 400 is fine for quick checks, and 1200+ for publication-quality plots.

---

## Output

### Console Output

```
══════════════════════════════════════════════════════════════
  CONTINUOUS BEAM ANALYSIS RESULTS
  Unit system: imperial
══════════════════════════════════════════════════════════════

── SUPPORT REACTIONS ─────────────────────────────────────
  x =     0.0000 in  (pinned)
      Vertical reaction:      38.7500 lb  (↑)
  x =    36.0000 in  (roller)
      Vertical reaction:     475.7812 lb  (↑)
  ...
```

For beams with elastic (spring) supports, the output also includes the displacement at each beam support:

```
  x =    60.0000 in  (beam)
      Vertical reaction:     115.3125 lb  (↑)
      Displacement:           -0.0231 in  (↓)
  ...

  Sum of vertical reactions: 800.0000 lb
  Total applied load:       800.0000 lb (downward)
  Equilibrium check:        1.180636e-05 lb

── EXTREME VALUES ───────────────────────────────────────
  Max shear force:         313.5937 lb
  Max bending moment:     2204.9765 lb·in
  Max deflection:          0.000369 in
  ...
```

The **equilibrium check** shows the difference between the sum of reactions and the total applied load. Values near zero (< 0.01) confirm the solution is correct.

### Diagram Output

The plot is saved to the specified path (default: `beam_analysis.png`) and contains:

1. **Beam Schematic** — free body diagram with:
   - Support symbols (triangle for pinned, triangle + circle for roller, hatched block for fixed, zigzag spring for beam supports)
   - Reaction forces with magnitude and direction (green ↑ for upward, orange ↓ for downward/uplift)
   - Applied loads with magnitudes and positions (red for point loads, blue for distributed, green for moments)

2. **Shear Force Diagram** — internal shear force V along the beam (in force units). Jumps occur at point loads and reactions.

3. **Bending Moment Diagram** — internal bending moment M along the beam (in force × length units). Positive = sagging (tension on bottom fiber), negative = hogging (tension on top fiber).

4. **Slope Diagram** — rotation angle θ of the beam cross-section (in radians). This is the first derivative of the deflection curve.

5. **Deflection Diagram** — vertical displacement δ of the beam (in length units). Negative = downward. Zero at pinned/roller/fixed support locations; non-zero at beam (spring) supports.

6. **Average Shear Stress Diagram** — τ_avg = V / A (in stress units). Same shape as the shear force diagram, scaled by cross-sectional area.

7. **Maximum Shear Stress Diagram** — τ_max = shear_shape_factor × V / A (in stress units). The actual peak shear stress at the neutral axis, accounting for cross-section shape. Compare this to the material's allowable shear stress.

8. **Bending Stress Diagram** — σ = M·c / I (in stress units). Maximum fiber stress at the extreme fiber of the cross-section. Compare this to the material's allowable bending stress.

---

## Programmatic Access

The solver can also be used as a library in your own scripts:

```python
from beam_solver import (
    ContinuousBeamSolver, BeamSection, Support,
    PointLoad, DistributedLoad, AppliedMoment
)

solver = ContinuousBeamSolver(
    length=120.0,
    section=BeamSection(E=29e6, I=30.8, A=2.96, c=3.94),
    supports=[
        Support(0, "pinned"),
        Support(60, "beam", stiffness=5000.0),   # elastic spring support
        Support(120, "roller"),
    ],
    distributed_loads=[
        DistributedLoad(20, 100, 10.0, 10.0),
    ],
)

solver.solve()

# Access results directly
print(solver.reactions)          # dict of {position: {vertical, moment, displacement*}}
print(solver.V_plot)             # shear force array
print(solver.M_plot)             # bending moment array
print(solver.delta_plot)         # deflection array
print(solver.theta_plot)         # slope array
print(solver.tau_avg_plot)       # average shear stress array
print(solver.tau_max_plot)       # maximum shear stress array (shape-factor adjusted)
print(solver.sigma_plot)         # bending stress array
print(solver.x_plot)             # x-coordinate array for all plots
# * 'displacement' key is only present for beam (spring) supports

solver.plot(save_path="my_beam.png")   # save plot to file
solver.plot()                          # display interactively
```

---

## Sign Conventions

| Quantity | Positive Direction |
|----------|-------------------|
| Point load magnitude | Downward ↓ |
| Distributed load intensity | Downward ↓ |
| Applied moment magnitude | Counterclockwise ↺ |
| Reaction (vertical) | Upward ↑ |
| Reaction (beam support) | Upward ↑ (R = k × δ, opposing displacement) |
| Reaction (moment, fixed) | Counterclockwise ↺ |
| Bending moment | Sagging (tension on bottom) |
| Deflection | Upward ↑ |
| Shear force | Standard beam convention (positive = clockwise couple on element) |

A **negative vertical reaction** indicates uplift — that support must be anchored to prevent the beam from lifting off.

---

## Theory

The solver uses the **direct stiffness method** with Euler-Bernoulli beam elements. Each element has two nodes with two degrees of freedom each (vertical displacement and rotation), connected by cubic Hermite shape functions.

**Assembly and solution:**
1. The beam is meshed into elements with nodes placed at all support, load, and endpoint locations.
2. Element stiffness matrices (4×4) are assembled into a global system.
3. Equivalent nodal forces are computed for distributed loads.
4. Boundary conditions are applied by constraining DOFs at support locations. For beam (spring) supports, the spring stiffness is added to the diagonal of the global stiffness matrix instead of constraining the DOF, so the displacement is solved for rather than prescribed as zero.
5. The reduced system Kf · Uf = Ff is solved for unknown displacements.
6. Reactions are recovered from R = K · U - F for rigid supports. For beam supports, the reaction is R = -k · δ (spring force opposing the displacement).

**Post-processing:**
- Deflection and slope are read directly from the solution vector.
- Bending moment is computed from the second derivative of the displacement field: M = EI · v″, evaluated using Hermite shape function derivatives.
- Shear force is computed from the third derivative: V = EI · v‴ (constant per element).
- Stresses are computed from the internal forces: σ = M·c/I (max bending stress at extreme fiber), τ_avg = V/A (average shear stress), and τ_max = shear_shape_factor × V/A (maximum shear stress at neutral axis).

---

## Limitations

- **Euler-Bernoulli theory** — assumes plane sections remain plane and neglects shear deformation. Accurate for slender beams (length/depth > 10). For deep beams, Timoshenko beam theory would be more appropriate.
- **Prismatic beam** — the cross-section (E, I, A, c) is constant along the entire length. Stepped or tapered beams are not supported.
- **Linear elastic** — assumes small deflections and linear material behavior. Not suitable for plastic analysis or large-deformation problems.
- **Shear stress approximation** — the maximum shear stress uses a single shape factor applied to V/A. For complex cross-sections (e.g. aluminum extrusions, multi-cell profiles), a more detailed analysis (VQ/It) may be needed. The shape factor approach is conservative for standard shapes.
- **2D analysis** — vertical loads and in-plane bending only. No lateral-torsional buckling, biaxial bending, or axial loads.

---

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| "Stiffness matrix is singular" | Beam is a mechanism (insufficient supports) | Add supports to prevent rigid-body motion |
| Equilibrium check > 0.1 | Mesh too coarse or numerical issues | Increase `n_elements` |
| Reactions don't match hand calculations | Check sign conventions | Positive loads = downward; positive moments = CCW |
| Plot looks jagged | Mesh too coarse | Increase `n_elements` to 800+ |
| Very large deflections | Units mismatch | Verify E, I, and length are in consistent units |
| "requires a positive stiffness value" | Beam support missing stiffness | Add `stiffness=<value>` to the `Support(kind="beam", ...)` |
| Beam support deflection too large | Spring stiffness too low | Increase the `stiffness` value or check units (lb/in vs N/mm) |

---

## License

This tool is provided as-is for educational and engineering use. The user is responsible for verifying results and ensuring they are appropriate for their application. This tool is not a substitute for professional structural engineering judgment.
