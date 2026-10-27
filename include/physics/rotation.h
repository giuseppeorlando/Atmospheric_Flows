/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022-2026 Giuseppe Orlando
 *
 * This code is free software; you can use it, redistribute it,
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * ------------------------------------------------------------------------
 *
 * Authors: Letizia Bottani, Giuseppe Orlando, 2026
 */
#pragma once

// @sect{Include files}

// We start by including the necessary header files
//
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

#include "../equation_data.h"

#include <string>

// @sect{Rotation}

// In this namespace, we implement the auxiliary structures and functions to take into account
// rotational effects through the Coriolis force on an f-plane. The rotation axis is aligned
// with the vertical direction, so that the Coriolis force acts only on the horizontal velocity components
// (see include/physics/vertical_slice.h for the relation between velocity components and directions of the mesh).
// The different formulations are discussed in https://arxiv.org/abs/2605.29598
//
namespace Rotation {
  using namespace dealii;

  /**
   * Available formulations for the rotation term:
   * - none: no rotation;
   * - R1: the Coriolis force at the current stage is treated through the fixed point loop;
   * - R2: the Coriolis force at the current stage is included in the linear system for the pressure
   *       through (A + R)^{-1} = (I + beta J)^{-1} A^{-1}, with beta = a_tilde*dt/Ro (see inverse_rotation_term).
   *       Moreover, also the gravity term of the energy equation at the current stage is included
   *       in the linear system for the pressure.
   *       WARNING: possible inconsistency in presence of a non-zero background velocity u_bar.
   *       The Coriolis force acts on the velocity relative to the background velocity, so that,
   *       at the current stage, the momentum equation reads
   *         rho*u + a_tilde*dt*(1/Ro)*rho*(k x (u - u_bar)) = rhs.
   *       Since u_bar is known, only the part acting on u has to be treated implicitly through
   *       (I + beta J), whereas the known contribution a_tilde*dt*(1/Ro)*rho*(k x u_bar) has to be
   *       moved to the rhs:
   *         rho*u + a_tilde*dt*(1/Ro)*rho*(k x u) = rhs + a_tilde*dt*(1/Ro)*rho*(k x u_bar).
   *       This contribution is currently not included (as in the original implementation), so that
   *       the Coriolis force at the current stage acts on u instead of u - u_bar. Hence, the formulation
   *       is consistent only for a background state at rest (e.g. inertia-gravity waves), whereas for
   *       u_bar != 0 (e.g. flows over mountains) a spurious Coriolis force is introduced.
   *       Possible solution: add, for the R2 formulation, the known contribution at the current stage
   *       to the rhs of the momentum equation (in EULEROperator::assemble_rhs_cell_term_momentum,
   *       next to the corresponding contribution of the R1 formulation):
   *         if(rotation == Rotation::Formulation::R2) {
   *           rotation_term -= a_tilde[IMEX_stage - 1][IMEX_stage - 1]*dt*
   *                            Rotation::coriolis_term(rho_s, phi_u_bar.get_value(q), inv_Ro);
   *         }
   *       ('rotation_term' is subtracted from the rhs, hence the minus sign). With this correction,
   *       R1 and R2 should coincide at convergence of the fixed point loop also for u_bar != 0.
   */
  enum class Formulation {none, R1, R2};

  /**
   * Convert the name of the formulation (as read from the parameter file) into the corresponding enumeration
   * @param name name of the formulation
   * @return the corresponding formulation
   */
  inline Formulation string_to_formulation(const std::string& name) {
    if(name == "none") {
      return Formulation::none;
    }
    else if(name == "R1") {
      return Formulation::R1;
    }
    else if(name == "R2") {
      return Formulation::R2;
    }
    else {
      AssertThrow(false, ExcMessage("Unknown formulation for the rotation term: " + name));
    }

    return Formulation::none;
  }

  /**
   * Coriolis force (non-dimensional) on an f-plane, i.e. (1/Ro)*rho*(k x u_rel), with k the unit vector
   * along the vertical direction. The Coriolis force acts on the velocity relative to the background velocity
   * @param rho density
   * @param u_rel velocity relative to the background velocity
   * @param inv_Ro inverse of the Rossby number
   * @return the Coriolis force
   */
  template<typename Number, typename T>
  inline DEAL_II_ALWAYS_INLINE
  Tensor<1, EquationData::n_comp_u, Number> coriolis_term(const Number& rho,
                                                          const Tensor<1, EquationData::n_comp_u, Number>& u_rel,
                                                          const T inv_Ro) {
    Tensor<1, EquationData::n_comp_u, Number> res;
    res[0] = -inv_Ro*rho*u_rel[1];
    res[1] = inv_Ro*rho*u_rel[0];

    return res;
  }

  /**
   * Inverse of (I + beta J), with J the matrix associated to k x u, i.e. J(u, v, w) = (-v, u, 0).
   * Since J acts only among the velocity components, this can be applied independently for each degree of freedom
   * @param v velocity components associated to a degree of freedom
   * @param beta coefficient of J (a_tilde*dt/Ro)
   * @return (I + beta J)^{-1} v
   */
  template<typename Number, typename T>
  inline DEAL_II_ALWAYS_INLINE
  Tensor<1, EquationData::n_comp_u, Number> inverse_rotation_term(const Tensor<1, EquationData::n_comp_u, Number>& v,
                                                                  const T beta) {
    const auto s = static_cast<T>(1.0)/(static_cast<T>(1.0) + beta*beta);

    Tensor<1, EquationData::n_comp_u, Number> res;
    res[0] = s*(v[0] + beta*v[1]);
    res[1] = s*(v[1] - beta*v[0]);
    res[2] = v[2];

    return res;
  }
} // namespace Rotation
