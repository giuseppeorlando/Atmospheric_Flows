# Atmospheric Flow Solver with Rotation

This repository is based on the atmospheric flow solver originally developed by
[Giuseppe Orlando](https://github.com/giuseppeorlando/Atmospheric_Flows/tree/Mesoscale).

The present branch extends the original framework to include rotational effects
in the simulation of atmospheric flows. The solver is based on a high-order
Discontinuous Galerkin formulation and is designed for the numerical simulation
of atmospheric flow problems.

The repository contains the source code, parameter files, and numerical test
cases associated with the present work. The parameter files required to run
the different test cases are provided in the `examples/` directory.

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

---

### Initial and Background Conditions

The initial atmospheric state is defined by

- Buoyancy frequency: `N`
- Background flow and thermodynamic conditions:
  - Background velocity: `u_bar`
  - Pressure at ground level (`z = 0`): `p_bar`
  - Temperature at ground level (`z = 0`): `T_bar`

---

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

---

## Non-Dimensional Formulation

The solver supports a fully non-dimensional formulation through the
specification of reference quantities:

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

## Running the Solver

Two main solver implementations are provided in the repository:

- `Atmospheric_Flow.cc` implements the three-dimensional solver, including
  simulations over orography.
- `Atmospheric_Flow_IGW.cc` implements the two-dimensional vertical-slice
  solver without orography. It is used, for example, for the
  two-dimensional inertia-gravity wave tests.

The corresponding executables can be compiled using the provided
`CMakeLists.txt`.

After compilation, the 3D solver can be launched as:

```bash
./build/Atmospheric_Flow
```

If no input file is specified, the solver attempts to read:
```bash
parameter-file.prm
```
from the current directory.
A custom parameter file can be specified using the -p (or --param) option:
```bash
./build/Atmospheric_Flow -p <parameter_file>
```
For parallel execution using MPI:
```bash
mpirun -np <number_of_processes> ./build/Atmospheric_Flow \
    -p <parameter-file>
```

For example:
```bash
mpirun -np 20 ./build/Atmospheric_Flow \
    -p examples/Mountain_3D_smalldomain/parameter-file.prm
```
A minimal help message can be displayed with:
```bash
./build/Atmospheric_Flow -h
```

The different test cases can be run using the parameter files provided in the
examples/ directory.

## Rotation Operators
Two different operators including the rotation term are provided:
- euler_operator_rotation.h
- euler_operator_rotation_R2.h
The differences between the two formulations are discussed in https://arxiv.org/abs/2605.29598

The operator used by the solver can be selected by including the corresponding
header file in the main solver source code.

## Additional Modules
The subfolder LAMINAR_INVERSE_MASS_MATRICES_PERTURBATION contains an
alternative solver implementation based on a perturbation formulation
(as described, e.g., in
https://www.sciencedirect.com/science/article/abs/pii/S002199911730030X).