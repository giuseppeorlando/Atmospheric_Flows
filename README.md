# Atmospheric Flow Solver

Source code for atmospheric dynamics simulations in [*deal.II*](https://www.dealii.org/).
The solver is based on a high-order Discontinuous Galerkin (DG) spatial discretization combined with an IMEX time integration strategy, as described in

> G. Orlando, T. Benacchio, and L. Bonaventura. "An IMEX-DG solver for atmospheric dynamics
> simulations with adaptive mesh refinement". In: Journal of Computational and Applied Mathematics
> (2023), p. 115124.

The solver can simulate three-dimensional flows over idealized orography (**versiera of Agnesi** profile) and
two-dimensional vertical slices, optionally including rotational (Coriolis) effects. The rotation terms have been
developed by Letizia Bottani (see also the repository
[Atmospheric_Flows_Rotation](https://github.com/letiziabottani/Atmospheric_Flows_Rotation)).

The solver is highly configurable via parameter files, allowing flexible setup of physical, geometrical, and numerical aspects of the simulation.

---

## Requirements

- [*deal.II*](https://www.dealii.org/) (version 9.5.1 or later) configured with
  - `DEAL_II_WITH_MPI = ON`
  - `DEAL_II_WITH_P4EST = ON`
  - `DEAL_II_WITH_HDF5 = ON`
  - `DEAL_II_WITH_LAPACK = ON`

---

## Compilation

The executable is built with CMake and placed in the root of the repository:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDEAL_II_DIR=/path/to/deal.II
cmake --build build
```

The script `configure.sh` performs the same operations in `Debug` mode.

### Dimension of the problem

The dimension of the problem is fixed at compile time through the constant `dim` at the beginning of `main()`
in `Atmospheric_Flow.cc` (by default `dim = 3`). Two-dimensional configurations are interpreted as vertical slices
in the x-z plane: the velocity has always three components, the second one being the out-of-plane component
(needed, e.g., for the Coriolis force; see `include/physics/vertical_slice.h`).
The inertia-gravity waves test case requires `dim = 2`.

---

## Running the Solver

After compilation, the executable can be launched as

```bash
./Atmospheric_Flow
```

If no input file is specified, the solver attempts to read `parameter-file.prm` from the current directory.
A custom parameter file can be specified using the `-p` (or `--param`) option

```bash
./Atmospheric_Flow -p <parameter_file>
```

For parallel execution using MPI

```bash
mpirun -np <number_of_processes> ./Atmospheric_Flow -p <parameter_file>
```

A help message with the description of all the parameters can be displayed with

```bash
./Atmospheric_Flow -h
```

whereas the parameters of the test case selected in the parameter file can be displayed with

```bash
./Atmospheric_Flow -p <parameter_file> --help_test_case
```

---

## Test cases

The test case is selected through the parameter `test_case` (subsection `Test case`), whereas its specific parameters
are read from the file indicated by `test_case_parameters`.

| `test_case` | Description | Dimension | Implementation |
|---|---|---|---|
| `3D mountain` | Flow over a mountain (background state with constant buoyancy frequency or isothermal) | 3 | `include/test_cases/mountain/ic_3D_mountain.h` |
| `inertia-gravity waves` | Perturbation of the (potential) temperature (background state isothermal or with constant buoyancy frequency) | 2 | `include/test_cases/IGW/ic_IGW.h` |

The hydrostatic background states shared by the test cases are implemented in `include/test_cases/background_state.h`.

---

## Examples

The global parameter files of some configurations (`parameter-file_*.prm`) are provided in the folders of the test cases,
together with the parameter files of the test cases (`ic_*.prm`). The examples have to be run from the root of the repository, e.g.

```bash
mpirun -np 20 ./Atmospheric_Flow -p include/test_cases/mountain/parameter-file_3D_nonhydrostatic.prm
```

| Global parameter file | Test case | Dimension |
|---|---|---|
| `include/test_cases/IGW/parameter-file_IGW.prm` | Inertia-gravity waves | 2 |
| `include/test_cases/IGW/parameter-file_IGW_planetary_scale.prm` | Inertia-gravity waves at planetary scale | 2 |
| `include/test_cases/IGW/parameter-file_IGW_rotating.prm` | Inertia-gravity waves with rotation (constant buoyancy frequency) | 2 |
| `include/test_cases/mountain/parameter-file_3D_nonhydrostatic.prm` | Mountain, non-hydrostatic regime (small domain, constant buoyancy frequency) | 3 |
| `include/test_cases/mountain/parameter-file_3D_hydrostatic.prm` | Mountain, hydrostatic regime (larger domain, isothermal atmosphere) | 3 |

All these examples are dimensional simulations, for which the tolerance of the linear solvers is absolute
(see [Global parameters](#global-parameters)).

---

## Rotation

Rotational effects are taken into account through the Coriolis force on an f-plane, acting on the velocity
relative to the background velocity (see `include/physics/rotation.h`). They are controlled by the parameters
(subsection `Physical data`)

- `Rossby`: the Rossby number;
- `rotation_formulation`: the formulation of the rotation term
  - `none`: no rotation (default);
  - `R1`: the Coriolis force at the current stage is treated through the fixed point loop;
  - `R2`: the Coriolis force at the current stage is included in the linear system for the pressure,
    whereas the known contribution of the background velocity is moved to the right-hand side.

At convergence of the fixed point loop, R1 and R2 solve the same discrete system.
The two formulations are discussed in https://arxiv.org/abs/2605.29598.

> **Note:** for a non-zero background velocity, the work of the Coriolis force acting on the velocity relative
> to the background velocity is included in the energy equation (see `include/physics/rotation.h`).

---

## Configuration

### Global parameters

The global parameters (file passed with `-p`) include

- Domain bounds:
  - `x_min`, `x_max`
  - `y_min`, `y_max` (only in three dimensions)
  - `z_min`, `z_max`
- Non-dimensional numbers: `Mach`, `Froude`, `Rossby`
- Reference quantities: `L_ref`, `u_ref`, `p_ref`, `rho_ref` (see [Non-Dimensional Formulation](#non-dimensional-formulation))
- Time step `dt` (or Courant number `CFL`), initial and final time
- Initial mesh resolution:
  - `n_elements_x`
  - `n_elements_y` (only in three dimensions)
  - `n_elements_z`
  - `n_of_refines` (number of global refinements)
- Solver tolerances:
  - Fixed-point iteration convergence: `||dp|| < atol_fixed_point + rtol_fixed_point * ||p||`
  - Linear solver convergence: `||Ax - b|| < atol_iterative + rtol_iterative * ||b||`.
    For dimensional simulations, `||b||` is dominated by the hydrostatic background state, so that
    `rtol_iterative` has to be very small and `atol_iterative` determines the accuracy (see the examples)

### Mountain test cases

- Orography parameters:
  - Mountain height: `h`
  - Mountain center coordinates: `xc`, `yc`
  - Mountain semi-width: `ac`
- Initial and background conditions:
  - Type of background state: `background` (`constant N` or `isothermal`)
  - Buoyancy frequency: `N` (only for `constant N`, `N = 0` for a neutral atmosphere)
  - Background velocity: `u_bar`
  - Pressure at ground level (`z = 0`): `p_bar`
  - Temperature (`isothermal`) or potential temperature at ground level (`constant N`): `T_bar`

The hydrostatic or non-hydrostatic regime is determined by the semi-width of the mountain with respect to `u_bar/N`,
rather than by the background state.
- Rayleigh damping layers (subsection `Rayleigh damping data`), configured independently along each spatial direction:
  - Vertical damping: `z_start`, `lambda_z`
  - Horizontal damping (x-direction): `x_start_left`, `lambda_x_left`, `x_start_right`, `lambda_x_right`
  - Horizontal damping (y-direction): `y_start_left`, `lambda_y_left`, `y_start_right`, `lambda_y_right`

### Inertia-gravity waves test case

- Perturbation of the temperature (`isothermal`) or of the potential temperature (`constant N`):
  - Amplitude: `dT` (`dT = 0` for a background state at rest)
  - Vertical scale (typically the height of the domain): `H`
  - Center: `xc`
  - Semi-width: `ac`
- Background state:
  - Type of background state: `background` (`isothermal` or `constant N`)
  - Buoyancy frequency: `N` (only for `constant N`, `N = 0` for a neutral atmosphere)
  - Background velocity: `u_bar`
  - Pressure at ground level (`z = 0`): `p_bar`
  - Temperature (`isothermal`) or potential temperature at ground level (`constant N`): `T_bar`

---

## Non-Dimensional Formulation

The solver supports a fully non-dimensional formulation through the specification of reference quantities

- Reference length: `L_ref`
- Reference velocity: `u_ref`
- Reference pressure: `p_ref`
- Reference density: `rho_ref`

> **Note:**
> - The **Mach**, **Froude** and **Rossby** numbers must still be provided explicitly and are enforced during the simulation.
> - These reference quantities do not necessarily correspond to physical scales; they are used purely for numerical normalization.
> - A dimensional simulation is obtained by setting all the reference quantities and the Mach number equal to 1
>   and the Froude number equal to 1/sqrt(g).

---

## Output

The solution (density `rho`, velocity `u`, pressure `p`) is saved in HDF5/XDMF and VTU format, together with the
background state and the potential temperature `theta`. The potential temperature is computed in Kelvin from the
pressure and the density expressed in physical units through the reference quantities, so that it is meaningful
both for dimensional and non-dimensional simulations. For two-dimensional configurations, the velocity components are
saved separately as `u`, `v` (out-of-plane), and `w`.

---

## Additional Modules

- `Atmospheric_Flow_turbulent.cc` and `turbulent_operator.h`: solver including turbulent vertical diffusion.
- `LAMINAR_INVERSE_MASS_MATRICES`, `LAMINAR_ITERATIVE_SOLVERS`: alternative implementations of the laminar solver.
- `LAMINAR_INVERSE_MASS_MATRICES_PERTURBATION`: alternative solver implementation based on a perturbation formulation
  (as described, e.g., in https://www.sciencedirect.com/science/article/abs/pii/S002199911730030X).
- `TURBULENT_INVERSE_MASS_MATRICES`, `TURBULENT_ITERATIVE_SOLVERS`: alternative implementations of the turbulent solver.

These modules do not include the rotation term.

---

## License

This code is distributed under the GNU Lesser General Public License v2.1
(see `LICENSE.txt`).
