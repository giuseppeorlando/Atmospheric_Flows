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
#include "../equation_data.h"

#include <cmath>

// @sect{Background states}

// In this namespace, we implement the vertical profiles of the hydrostatic background states
// shared by several test cases:
// - 'isothermal': isothermal atmosphere (constant temperature T_bar);
// - 'constant N': atmosphere with constant buoyancy frequency N and potential temperature T_bar at z = 0
//   (neutral atmosphere for N = 0).
// The Exner pressure is normalized with its value at z = 0, so that the pressure reads p_bar*pi^(1/Gamma),
// with p_bar the pressure at z = 0 and Gamma = (gamma - 1)/gamma. The height is z = x_vertical*L_ref,
// with x_vertical the (possibly non-dimensional) vertical coordinate.
//
namespace BackgroundState {
  /**
   * Exner pressure (normalized with its value at z = 0) of an isothermal atmosphere
   * @param x_vertical vertical coordinate
   * @param L_ref reference length
   * @param T_bar background temperature
   * @return the Exner pressure
   */
  template<typename T>
  T exner_isothermal(const T x_vertical, const T L_ref, const T T_bar) {
    const auto Gamma = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                       static_cast<T>(EquationData::Cp_Cv);
    const auto Cp    = static_cast<T>(EquationData::R)/Gamma;

    return std::exp(-static_cast<T>(EquationData::g)/(Cp*T_bar)*x_vertical*L_ref);
  }

  /**
   * Exner pressure (normalized with its value at z = 0) of an atmosphere with constant buoyancy frequency
   * @param x_vertical vertical coordinate
   * @param L_ref reference length
   * @param N buoyancy frequency (N = 0 for a neutral atmosphere)
   * @param T_bar potential temperature at z = 0
   * @return the Exner pressure
   */
  template<typename T>
  T exner_constant_N(const T x_vertical, const T L_ref, const T N, const T T_bar) {
    const auto Gamma = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                       static_cast<T>(EquationData::Cp_Cv);

    // Neutral atmosphere (limit for N -> 0)
    if(N == static_cast<T>(0.0)) {
      return static_cast<T>(1.0)
           - static_cast<T>(EquationData::g)*Gamma/(static_cast<T>(EquationData::R)*T_bar)*x_vertical*L_ref;
    }

    return static_cast<T>(1.0)
         - static_cast<T>(EquationData::g)*static_cast<T>(EquationData::g)/(N*N)*
           Gamma/(static_cast<T>(EquationData::R)*T_bar)*
           (static_cast<T>(1.0) - std::exp(-N*N/static_cast<T>(EquationData::g)*x_vertical*L_ref));
  }

  /**
   * Potential temperature of an atmosphere with constant buoyancy frequency
   * @param x_vertical vertical coordinate
   * @param L_ref reference length
   * @param N buoyancy frequency (N = 0 for a neutral atmosphere)
   * @param T_bar potential temperature at z = 0
   * @return the potential temperature
   */
  template<typename T>
  T potential_temperature_constant_N(const T x_vertical, const T L_ref, const T N, const T T_bar) {
    return T_bar*std::exp(N*N/static_cast<T>(EquationData::g)*x_vertical*L_ref);
  }
} // namespace BackgroundState
