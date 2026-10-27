# Atmospheric Flow Solver with Rotation

Source code for atmospheric dynamics simulations in [*deal.II*](https://www.dealii.org/).
The solver is based on a high-order Discontinuous Galerkin (DG) spatial discretization combined with an IMEX time integration strategy, as described in

> G. Orlando, T. Benacchio, and L. Bonaventura. "An IMEX-DG solver for atmospheric dynamics
> simulations with adaptive mesh refinement". In: Journal of Computational and Applied Mathematics
> (2023), p. 115124.

The present version extends the solver over orography to include rotational (Coriolis) effects. The rotation terms have been developed by Letizia Bottani, starting from the solver developed by Giuseppe Orlando (branch `Mesoscale` of this repository). The rotational formulations are discussed in https://arxiv.org/abs/2605.29598.

---

## About this version

The code on the `master` branch corresponds to the branch `rotation` of this
repository. The same code is also available in the separate repository
[Atmospheric_Flows_Rotation](https://github.com/letiziabottani/Atmospheric_Flows_Rotation), which contains only the files related to the solver with rotation (it does not include the additional modules described below).

---

## Requirements

- [*deal.II*](https://www.dealii.org/) (version 9.5.1 or later) configured with
  - `DEAL_II_WITH_MPI = ON`
  - `DEAL_II_WITH_P4EST = ON`
  - `DEAL_II_WITH_HDF5 = ON`
  - `DEAL_II_WITH_LAPACK = ON`

---

## Repository structure

- `Atmospheric_Flow.cc`: three-dimensional solver, including simulations over orography.
- `Atmospheric_Flow_IGW.cc`: two-dimensional vertical-slice solver without orography, used, e.g., for the inertia-gravity wave tests.
- `euler_operator_rotation.h`, `euler_operator_rotation_R2.h`: Euler operators
  including the rotation term (see [Rotation operators](#rotation-operators)).
- `euler_operator.h`: Euler operator without rotation.
- `include/`: initial and boundary conditions, Rayleigh damping, mapping, numerical fluxes, runtime parameters, and time integrators.
- `examples/`: parameter files for the test cases.
- `parameter-file.prm`: default parameter file.

---

## Compilation

The executable is built with CMake. The source file to be compiled is selected
through the variable `TARGET` in `CMakeLists.txt` (by default `Atmospheric_Flow`):

```bash
mkdir build && cd build
cmake -DDEAL_II_DIR=/path/to/deal.II ..
make release
make
```

To build the two-dimensional vertical-slice solver, set `SET(TARGET "Atmospheric_Flow_IGW")` in `CMakeLists.txt`.

---

## Running the Solver

After compilation, the 3D solver can be launched as

```bash
./build/Atmospheric_Flow
```

If no input file is specified, the solver attempts to read
`parameter-file.prm` from the current directory.
A custom parameter file can be specified using the `-p` (or `--param`) option

```bash
./build/Atmospheric_Flow -p <parameter_file>
```

For parallel execution using MPI

```bash
mpirun -np <number_of_processes> ./build/Atmospheric_Flow -p <parameter_file>
```

For example

```bash
mpirun -np 20 ./build/Atmospheric_Flow -p examples/Mountain_3D_smalldomain/parameter-file.prm
```

A minimal help message can be displayed with

```bash
./build/Atmospheric_Flow -h
```

---

## Test cases

The parameter files for the test cases are provided in `examples/`:

| Folder | Description | Solver |
|---|---|---|
| `examples/IGW` | Inertia-gravity waves in a vertical slice | `Atmospheric_Flow_IGW` |
| `examples/IGW_planetary_scale` | Inertia-gravity waves in a vertical slice at planetary scale | `Atmospheric_Flow_IGW` |
| `examples/Mountain_3D_smalldomain` | Three-dimensional flow over orography (small domain) | `Atmospheric_Flow` |
| `examples/Mountain_3D_largerdomain` | Three-dimensional flow over orography (larger domain) | `Atmospheric_Flow` |

---

## Rotation operators

Two different operators including the rotation term are provided

- `euler_operator_rotation.h`
- `euler_operator_rotation_R2.h`

The differences between the two formulations are discussed in https://arxiv.org/abs/2605.29598.

The operator is selected in the main source file:

- `Atmospheric_Flow.cc` includes `euler_operator_rotation.h`;
- `Atmospheric_Flow_IGW.cc` selects the operator by defining exactly one of
  `ROTATION_R1` (`euler_operator_rotation.h`) or `ROTATION_R2`
  (`euler_operator_rotation_R2.h`) at the beginning of the file.

---

## Configuration

### Domain and Orography Parameters

The computational domain and terrain geometry, when applicable, can be defined
through

- Domain bounds:
  - `x_min`, `x_max`
  - `y_min`, `y_max`
  - `z_min`, `z_max`
- Orography parameters:
  - Mountain height: `h`
  - Mountain center coordinates: `xc`, `yc`
  - Mountain semi-width: `ac`

### Initial and Background Conditions

The initial atmospheric state is defined by

- Buoyancy frequency: `N`
- Background flow and thermodynamic conditions:
  - Background velocity: `u_bar`
  - Pressure at ground level (`z = 0`): `p_bar`
  - Temperature at ground level (`z = 0`): `T_bar`

### Numerical Parameters

#### Mesh

- Initial mesh resolution:
  - `n_elements_x`
  - `n_elements_y`
  - `n_elements_z`

#### Rayleigh Damping

Rayleigh damping layers can be configured independently along each spatial
direction:

- Vertical damping:
  - `z_start`, `lambda_z`
- Horizontal damping (x-direction):
  - `x_start_left`, `lambda_x_left`
  - `x_start_right`, `lambda_x_right`
- Horizontal damping (y-direction):
  - `y_start_left`, `lambda_y_left`
  - `y_start_right`, `lambda_y_right`

#### Solver Tolerances

- Fixed-point iteration convergence:
  - `||dp|| < atol_iterative + rtol_iterative * ||p||`
- Linear solver convergence:
  - `||Ax - b|| < atol_fixed_point + rtol_fixed_point * ||b||`

### Non-Dimensional Formulation

The solver supports a fully non-dimensional formulation through the specification of reference quantities:

- Reference length: `L_ref`
- Reference velocity: `u_ref`
- Reference pressure: `p_ref`
- Reference density: `rho_ref`

> **Note:**
> - The **Mach**, **Froude** and **Rossby** numbers must still be provided explicitly and
>   are enforced during the simulation.
> - These reference quantities do not necessarily correspond to physical
>   scales; they are used purely for numerical normalization.

---

## Additional Modules

The following modules are provided in separate folders and do not include the
rotation term:

- `LAMINAR_INVERSE_MASS_MATRICES`, `LAMINAR_ITERATIVE_SOLVERS`: alternative
  implementations of the laminar solver.
- `LAMINAR_INVERSE_MASS_MATRICES_PERTURBATION`: alternative solver implementation
  based on a perturbation formulation (as described, e.g., in
  https://www.sciencedirect.com/science/article/abs/pii/S002199911730030X).
- `TURBULENT_INVERSE_MASS_MATRICES`, `TURBULENT_ITERATIVE_SOLVERS`, together with
  `Atmospheric_Flow_turbulent.cc` and `turbulent_operator.h`: solver including
  turbulent vertical diffusion.

---

## License

This code is distributed under the GNU Lesser General Public License v2.1
(see `LICENSE.txt`).
