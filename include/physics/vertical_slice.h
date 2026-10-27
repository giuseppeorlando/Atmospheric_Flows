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
#include <deal.II/base/tensor.h>

#include "../equation_data.h"

#include <type_traits>

// @sect{Vertical slice}

// In this namespace, we implement auxiliary functions to deal with a velocity field
// with EquationData::n_comp_u components on a mesh of dimension dim.
// In three dimensions, the velocity components coincide with the directions of the mesh
// and all the functions reduce to the identity. In two dimensions, the mesh is interpreted
// as a vertical slice in the x-z plane: the velocity components 0 and 2 correspond to the
// two directions of the mesh, whereas the velocity component 1 is the out-of-plane one
// (needed, e.g., to take into account the Coriolis force).
//
namespace VerticalSlice {
  using namespace dealii;

  /**
   * Velocity component associated to a direction of the mesh
   * @param d direction of the mesh
   * @return index of the corresponding velocity component
   */
  template<unsigned dim>
  constexpr unsigned velocity_component(const unsigned d) {
    static_assert(dim == 2 || dim == EquationData::n_comp_u,
                  "Only two-dimensional (vertical slice) and three-dimensional configurations are supported");

    return (d == dim - 1) ? EquationData::n_comp_u - 1 : d;
  }

  /**
   * Type for tensors with one entry per velocity component and one per direction of the mesh
   * (e.g. the gradient of the velocity). It coincides with the type employed by FEEvaluation
   */
  template<unsigned dim, typename Number>
  using velocity_gradient_type = std::conditional_t<dim == EquationData::n_comp_u,
                                                    Tensor<2, dim, Number>,
                                                    Tensor<1, EquationData::n_comp_u, Tensor<1, dim, Number>>>;

  /**
   * Express a vector defined along the directions of the mesh (e.g. the unit normal vector)
   * in the space of velocity components
   * @param v vector along the directions of the mesh
   * @return the corresponding vector in the space of velocity components
   */
  template<unsigned dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE
  Tensor<1, EquationData::n_comp_u, Number> to_velocity_space(const Tensor<1, dim, Number>& v) {
    if constexpr(dim == EquationData::n_comp_u) {
      return v;
    }
    else {
      Tensor<1, EquationData::n_comp_u, Number> res;
      for(unsigned d = 0; d < dim; ++d) {
        res[velocity_component<dim>(d)] = v[d];
      }

      return res;
    }
  }

  /**
   * Restrict a vector in the space of velocity components to the directions of the mesh
   * @param u vector in the space of velocity components
   * @return the components of u along the directions of the mesh
   */
  template<unsigned dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE
  Tensor<1, dim, Number> to_mesh_space(const Tensor<1, EquationData::n_comp_u, Number>& u) {
    if constexpr(dim == EquationData::n_comp_u) {
      return u;
    }
    else {
      Tensor<1, dim, Number> res;
      for(unsigned d = 0; d < dim; ++d) {
        res[d] = u[velocity_component<dim>(d)];
      }

      return res;
    }
  }

  /**
   * Restrict the second index of a tensor in the space of velocity components to the directions of the mesh
   * (e.g. the momentum flux to be tested against the gradient of the test functions)
   * @param t tensor in the space of velocity components
   * @return the tensor with second index along the directions of the mesh
   */
  template<unsigned dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE
  velocity_gradient_type<dim, Number> to_mesh_space(const Tensor<2, EquationData::n_comp_u, Number>& t) {
    if constexpr(dim == EquationData::n_comp_u) {
      return t;
    }
    else {
      velocity_gradient_type<dim, Number> res;
      for(unsigned i = 0; i < EquationData::n_comp_u; ++i) {
        res[i] = to_mesh_space<dim>(t[i]);
      }

      return res;
    }
  }
} // namespace VerticalSlice
