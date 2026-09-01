# Lid-driven cavity

A finite-volume incompressible Navier–Stokes solver on a staggered grid, written from scratch in
MATLAB and validated against Ghia, Ghia and Shin (1982). The solver diverged at a time step that
both textbook stability conditions said was safe; finding out why is the part of this worth
reading.

Written for **ME 6434 Advanced CFD**, Virginia Tech, Spring 2022 (Prof. Danesh Tafti).

Fluid sits in a square cavity. Three walls are fixed; the top wall slides at constant velocity. The
question is how the u and v profiles develop as the flow reaches steady state.

![lid_driven_cavity](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/382fa46b-ac14-42aa-8618-fbe46c894d83)

The part worth reading is [Why the 32×32 grid diverged](#why-the-3232-grid-diverged), where the
solver blew up at a time step that both textbook stability conditions said was safe.

📄 [Technical report](Technical_report.pdf) · [Problem statement](Problem_statement.pdf) ·
[Ghia et al. (1982)](ghia1982.pdf) · [`Codes/lid_driven_cavity_2d.m`](Codes/lid_driven_cavity_2d.m)

## Formulation

Non-dimensionalised on the cavity length $L$ and lid velocity $u_T$, so $Re = u_T L / \nu$:

$$\nabla \cdot \mathbf{u} = 0$$

$$\frac{\partial \mathbf{u}}{\partial t} + \nabla\cdot(\mathbf{u}\mathbf{u}) = -\nabla p + \frac{1}{Re}\nabla^2\mathbf{u}$$

**Staggered grid.** Pressure at cell centres, velocities on cell faces. This falls out of the
continuity equation itself, writing $\nabla\cdot\mathbf{u}=0$ over a control volume needs the
velocities on the faces. A collocated arrangement admits the checkerboard pressure mode; staggering
removes it without needing Rhie–Chow interpolation. Momentum control volumes are therefore offset
from the scalar ones, staggered to the negative faces so velocity indices match the cell-centred
index.

**Discretisation.** 2nd-order central differences in space, 2nd-order Adams–Bashforth in time.

**Fractional step (Perot's block-LU form).** Integrating an incomplete momentum equation gives an
intermediate velocity that is not divergence-free; a pressure Poisson solve projects it onto the
divergence-free space without changing vorticity. Written as a block LU factorisation this
decouples into three steps per time step:

$$A\mathbf{u}^\ast = \mathbf{r}, \qquad DG\,p = D\mathbf{u}^\ast, \qquad \mathbf{u}^{n+1} = \mathbf{u}^\ast - Gp$$

## Implementation

| | |
|---|---|
| Grid | 128 × 128 (32² and 64² also run for the grid study) |
| Pressure Poisson | **SOR**, hand-written, residual tolerance 1e−5 |
| Poisson BCs | Neumann, imposed by zeroing the wall-side coefficients rather than using ghost values |
| Steady-state criterion | L2 norm of the velocity change per step < 1e−8 |
| Time step | not fixed: chosen from three stability limits, see below |

Everything is in one MATLAB file with no toolbox dependencies. The Poisson solve is written out
rather than handed to `\`, because the point of the exercise was the discretisation.

## Choosing the time step

Rather than hard-code Δt, the solver takes the most restrictive of three limits:

```matlab
dt_max = min([ dx/uT, ...                              % linear CFL
               2^(2/3)*C^(1/3)*(dx/uT)^(4/3), ...      % non-linear CFL (Adams-Bashforth)
               0.25*dx*dx/nu ]);                       % 2D Neumann (viscous)
```

The middle term is not standard, and the reason it is there is the interesting part.

## Why the 32×32 grid diverged

On the 32×32 grid the solver diverged at **Δt = 0.02**. It should not have:

| Condition | Limit at 32×32 |
|---|---|
| Linear CFL, $\Delta t \le \Delta x/u$ | 0.0313 |
| 2D Neumann, $\Delta t \le 0.25\,\Delta x^2/\nu$ | 0.0244 |

Δt = 0.02 is below both. The obvious conclusions (a coding error, or an unstable scheme) were
both wrong.

Higher-order explicit time integrators (Adams–Bashforth, Runge–Kutta) carry a **non-linear**
CFL-type restriction that is more stringent than the linear one for convection-dominated flow.
Schneider et al. showed this numerically; Deriaz derived it from a von Neumann analysis of the 2D
Burgers and Euler equations. For 2nd-order Adams–Bashforth:

$$\Delta t \le 2^{2/3} C^{1/3} \left(\frac{\Delta x}{u}\right)^{4/3}, \qquad 0 < C \le 1$$

At 32×32 with $C = 1$ this gives **Δt ≤ 0.0156**, below the 0.02 that failed. The divergence was
correct behaviour under a condition neither textbook criterion captures.

That derived limit is the middle term in `dt_max` above, so the solver now picks a stable step on
its own at any resolution.

## Which limit binds depends on Re

Same 128×128 grid, two Reynolds numbers:

| Condition | Re = 100 | Re = 1000 |
|---|---|---|
| Linear CFL | 0.0078 | 0.0078 |
| 2D Neumann (viscous) | **0.0015** | 0.0153 |
| Non-linear CFL | 0.0025 | **0.0025** |

At Re = 100 the viscous limit binds; at Re = 1000 the non-linear CFL does. The viscous limit scales
as $\Delta x^2/\nu$ so it relaxes as Re rises, while the non-linear CFL is independent of Re. Which
is the physical statement that low-Re flow is stability-limited by diffusion and high-Re flow by
convection, so a solver hard-coded to the viscous limit would be needlessly slow at Re = 100 and
**unstable at Re = 1000**.

## Validation against Ghia et al. (1982)

At Re = 100, centreline velocities against the tabulated benchmark:

### u along the vertical centreline
![u velocity](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/f7b1d83b-5d50-4fd2-b818-5e9318253dd8)

### v along the horizontal centreline
![v velocity](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/867552d2-42df-46fe-9c01-5fb507c7a974)

The u profile matches closely. The v profile reproduces the oscillatory shape but the values agree
less well, reported as found rather than tuned away.

## Cost scaling

CPU time to convergence against grid resolution follows a power law $y = a x^{b}$, fitted in
MATLAB across 32², 64² and 128² at three time steps each. The exponent runs from **3.08** at the
smallest time step to **4.74** at the largest.

Above the ideal ~3 (cells × steps, with steps growing as the stability limit tightens), because the
SOR Poisson solve needs more sweeps per step as the grid refines: the iterative solve, not the
discretisation, is what makes refinement expensive here.

## Running it

Open [`Codes/lid_driven_cavity_2d.m`](Codes/lid_driven_cavity_2d.m), and run. Set `nx`, `ny` for
resolution and `nu` for Reynolds number: as committed it is 128 × 128 at `nu = 0.001`, i.e.
**Re = 1000**. The validation plots above are Re = 100, which is `nu = 0.01`.

No toolboxes required.

## References

1. K. Schneider, N. Kevlahan, M. Farge. *Comparison of an adaptive wavelet method and nonlinearly
   filtered pseudospectral methods for two-dimensional turbulence.* Theor. Comput. Fluid Dyn. 9
   (1997) 191–206.
2. E. Deriaz. *Stability conditions for the numerical solution of convection-dominated problems
   with skew-symmetric discretizations.* SIAM J. Numer. Anal. 50 (2012) 1058–1085.
3. U. Ghia, K. N. Ghia, C. T. Shin. *High-Re solutions for incompressible flow using the
   Navier–Stokes equations and a multigrid method.* J. Comput. Phys. 48 (1982) 387–411.
4. J. B. Perot. *An analysis of the fractional step method.* J. Comput. Phys. 108 (1993) 51–58.
