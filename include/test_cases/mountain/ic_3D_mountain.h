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
#include "../background_state.h"
#include "../test_case.h"

#include "../../equation_data.h"

#include <cmath>

// @sect{Initial conditions}

// In this namespace, we declare the initial background conditions for a flow over a mountain
// with a uniform horizontal background velocity. Two background states are available
// (see include/test_cases/background_state.h):
// - 'constant N': atmosphere with constant buoyancy frequency (neutral for N = 0);
// - 'isothermal': isothermal atmosphere.
// Notice that the hydrostatic or non-hydrostatic regime is determined by the semi-width
// of the mountain with respect to u_bar/N, rather than by the background state.
// Since the initial conditions coincide with the background state, the test case does not provide
// the background state separately (see TestCaseBase in include/test_cases/test_case.h).
// The last coordinate is the vertical one
//
namespace IC_Mountain3D {
  using namespace dealii;

  /**
   * We declare now the class that describes the initial condition for the velocity.
   */
  template<unsigned dim, typename T = double>
  class Velocity: public Function<dim, T> {
  public:
    /**
     * Class constructor
     * @param parameters_ auxiliary structure with parameters specific of the test case
     * @param data_ auxiliary structure with 'generic' parameters
     */
    Velocity(RunTimeParameters::MountainParameters& parameters_,
             RunTimeParameters::Data_Storage& data_);

    /**
     * Evaluation of the velocity for each component
     * @param p point coordinate
     * @param component spatial component to be evaluated
     */
    virtual T value(const Point<dim, T>& p,
                    const unsigned       component = 0) const override;

    /**
     * Vector evaluation of the velocity
     * @param p point coordinate
     * @return values vector with all components of the velocity
     */
    virtual void vector_value(const Point<dim, T>& p,
                              Vector<T>&           values) const override;

  private:
    RunTimeParameters::MountainParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;             /*!< 'Global' parameters */
  };

  // Constructor which relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Velocity<dim, T>::Velocity(RunTimeParameters::MountainParameters& parameters_,
                             RunTimeParameters::Data_Storage& data_):
    Function<dim, T>(EquationData::n_comp_u, data_.initial_time), parameters(parameters_), data(data_) {}

  // Specify the value for each spatial component. This function is overriden.
  //
  template<unsigned dim, typename T>
  T Velocity<dim, T>::value(const Point<dim, T>& p,
                            const unsigned       component) const {
    AssertIndexRange(component, EquationData::n_comp_u);

    if(component == 0) {
      return parameters.u_bar/data.u_ref;
    }
    else {
      return static_cast<T>(0.0);
    }
  }

  // Put together for a vector evalutation of the velocity.
  //
  template<unsigned dim, typename T>
  void Velocity<dim, T>::vector_value(const Point<dim, T>& p,
                                      Vector<T>&           values) const {
    Assert(values.size() == EquationData::n_comp_u, ExcDimensionMismatch(values.size(), EquationData::n_comp_u));

    for(unsigned i = 0; i < EquationData::n_comp_u; ++i) {
      values[i] = value(p, i);
    }
  }


  /**
   * We do the same for the pressure.
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

  // Constructor which again relies on the 'Function' constructor.
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

    const auto pi_bar = (parameters.background == "isothermal") ?
                        BackgroundState::exner_isothermal<T>(p[dim - 1], data.L_ref, parameters.T_bar) :
                        BackgroundState::exner_constant_N<T>(p[dim - 1], data.L_ref, parameters.N, parameters.T_bar);

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

    // Background density (at z = 0)
    const auto rho_bar = parameters.p_bar/(static_cast<T>(EquationData::R)*parameters.T_bar);

    // Isothermal atmosphere
    if(parameters.background == "isothermal") {
      const auto Gamma  = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                          static_cast<T>(EquationData::Cp_Cv);

      const auto pi_bar = BackgroundState::exner_isothermal<T>(p[dim - 1], data.L_ref, parameters.T_bar);

      return (rho_bar/data.rho_ref)*std::pow(pi_bar, static_cast<T>(1.0)/Gamma);
    }

    // Atmosphere with constant buoyancy frequency
    const auto pi_bar    = BackgroundState::exner_constant_N<T>(p[dim - 1], data.L_ref, parameters.N, parameters.T_bar);
    const auto theta_bar = BackgroundState::potential_temperature_constant_N<T>(p[dim - 1], data.L_ref,
                                                                                parameters.N, parameters.T_bar);

    return (rho_bar/data.rho_ref)*
           parameters.T_bar/theta_bar*std::pow(pi_bar, static_cast<T>(1.0)/
                                               (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0)));
  }

} // namespace IC_Mountain3D

/**
 * @brief 3D mountain test case
 *
 * This type owns all mountain-specific configuration and hides the concrete
 * implementations from the solver.
 */
template<unsigned dim, typename T = double>
using Mountain3DTestCase = TestCase<dim, T,
                                    RunTimeParameters::MountainParameters,
                                    IC_Mountain3D::Density<dim, T>,
                                    IC_Mountain3D::Velocity<dim, T>,
                                    IC_Mountain3D::Pressure<dim, T>>;
