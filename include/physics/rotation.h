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
   * - R1: the Coriolis force at the current stage is treated through the fixed point loop.
   */
  enum class Formulation {none, R1};

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
} // namespace Rotation
