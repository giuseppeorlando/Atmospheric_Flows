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

// We start by including the necessary deal.II header files and some C++
// related ones
//
#include <deal.II/base/function.h>

#include "mountain_parameters.h"
#include "ic_3D_nonhydrostatic.h"
#include "../test_case.h"

#include "../../equation_data.h"

#include <cmath>

// @sect{Initial conditions}

// In this namespace, we declare the initial background conditions for a flow over a mountain
// in an isothermal atmosphere (typically employed for hydrostatic configurations).
// The initial velocity is the same of the non-hydrostatic test case
// (see ic_3D_nonhydrostatic.h), whereas the buoyancy frequency is not employed,
// since it is determined by the background temperature
//
namespace IC_Hydrostatic3D {
  using namespace dealii;

  /**
   * We declare now the class that describes the initial condition for the pressure.
   */
  template<unsigned dim, typename T = double>
  class Pressure: public Function<dim, T> {
  public:
    /**
     * Class constructor
     * @param parameters_ auxiliary structure with parameters specific of the test case
     * @param data_ auxiliary structure with 'generic' parameters
     */
    Pressure(RunTimeParameters::MountainParameters& parameters_,
             RunTimeParameters::Data_Storage& data_);

    /**
     * Evaluation of the pressure
     * @param p point coordinate
     * @param component spatial component to be evaluated (unused)
     */
    virtual T value(const Point<dim, T>& p,
                    const unsigned       component = 0) const override;

  private:
    RunTimeParameters::MountainParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;             /*!< 'Global' parameters */
  };

  // Constructor which relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Pressure<dim, T>::Pressure(RunTimeParameters::MountainParameters& parameters_,
                             RunTimeParameters::Data_Storage& data_):
    Function<dim, T>(1, data_.initial_time), parameters(parameters_), data(data_) {}

  // Evaluation depending on the spatial coordinates. The input argument 'component'
  // will be unused but it has to be kept to override
  //
  template<unsigned dim, typename T>
  T Pressure<dim, T>::value(const Point<dim, T>& p,
                            const unsigned       component) const {
    (void)component;
    AssertIndexRange(component, 1);

    const auto Gamma  = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                        static_cast<T>(EquationData::Cp_Cv);
    const auto Cp     = static_cast<T>(EquationData::R)/Gamma;

    const auto pi_bar = std::exp(-static_cast<T>(EquationData::g)/(Cp*parameters.T_bar)*p[dim - 1]*data.L_ref);

    return (parameters.p_bar/data.p_ref)*std::pow(pi_bar, static_cast<T>(1.0)/Gamma);
  }


  /**
   * We do the same for the density.
   */
  template<unsigned dim, typename T = double>
  class Density: public Function<dim, T> {
  public:
    /**
     * Class constructor
     * @param parameters_ auxiliary structure with parameters specific of the test case
     * @param data_ auxiliary structure with 'generic' parameters
     */
    Density(RunTimeParameters::MountainParameters& parameters_,
            RunTimeParameters::Data_Storage& data_);

    /**
     * Evaluation of the density
     * @param p point coordinate
     * @param component spatial component to be evaluated (unused)
     */
    virtual T value(const Point<dim, T>& p,
                    const unsigned       component = 0) const override;

  private:
    RunTimeParameters::MountainParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;             /*!< 'Global' parameters */
  };

  // Constructor which again relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Density<dim, T>::Density(RunTimeParameters::MountainParameters& parameters_,
                           RunTimeParameters::Data_Storage& data_):
    Function<dim, T>(1, data_.initial_time), parameters(parameters_), data(data_) {}

  // Evaluation depending on the spatial coordinates. The input argument 'component'
  // will be unused but it has to be kept to override
  //
  template<unsigned dim, typename T>
  T Density<dim, T>::value(const Point<dim, T>& p,
                           const unsigned       component) const {
    (void)component;
    AssertIndexRange(component, 1);

    const auto Gamma  = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                        static_cast<T>(EquationData::Cp_Cv);
    const auto Cp     = static_cast<T>(EquationData::R)/Gamma;

    const auto pi_bar = std::exp(-static_cast<T>(EquationData::g)/(Cp*parameters.T_bar)*p[dim - 1]*data.L_ref);

    const auto rho_bar = parameters.p_bar/(static_cast<T>(EquationData::R)*parameters.T_bar);

    return (rho_bar/data.rho_ref)*std::pow(pi_bar, static_cast<T>(1.0)/Gamma);
  }

} // namespace IC_Hydrostatic3D

/**
 * @brief 3D hydrostatic test case
 *
 * This type owns all 3D hydrostatic-specific configuration and hides the concrete
 * implementations from the solver.
 */
template<unsigned dim, typename T = double>
using Hydrostatic3DTestCase = TestCase<dim, T,
                                       RunTimeParameters::MountainParameters,
                                       IC_Hydrostatic3D::Density<dim, T>,
                                       IC_NonHydrostatic3D::Velocity<dim, T>,
                                       IC_Hydrostatic3D::Pressure<dim, T>>;
