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
#include <deal.II/base/numbers.h>

#include "IGW_parameters.h"
#include "../background_state.h"
#include "../test_case.h"

#include "../../equation_data.h"

#include <cmath>

// @sect{Initial conditions}

// In this namespace, we declare the initial conditions for the inertia-gravity waves
// test case: a perturbation superimposed to a background state at rest (or with a uniform
// horizontal background velocity). Two background states are available:
// - 'isothermal': isothermal atmosphere with a perturbation of the temperature;
// - 'constant N': atmosphere with constant buoyancy frequency (neutral for N = 0) with a perturbation
//   of the potential temperature (see include/test_cases/background_state.h).
// For a zero amplitude of the perturbation, the background state at rest is obtained
// (e.g. to test the preservation of hydrostatic equilibrium). The last coordinate
// is the vertical one (see include/physics/vertical_slice.h)
//
namespace IC_IGW {
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
    Velocity(RunTimeParameters::IGWParameters& parameters_,
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
    RunTimeParameters::IGWParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;        /*!< 'Global' parameters */
  };

  // Constructor which relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Velocity<dim, T>::Velocity(RunTimeParameters::IGWParameters& parameters_,
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
    Pressure(RunTimeParameters::IGWParameters& parameters_,
             RunTimeParameters::Data_Storage& data_);

    /**
     * Evaluation of the pressure
     * @param p point coordinate
     * @param component spatial component to be evaluated (unused)
     */
    virtual T value(const Point<dim, T>& p,
                    const unsigned       component = 0) const override;

  private:
    RunTimeParameters::IGWParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;        /*!< 'Global' parameters */
  };

  // Constructor which again relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Pressure<dim, T>::Pressure(RunTimeParameters::IGWParameters& parameters_,
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

    // Background state with constant buoyancy frequency
    if(parameters.background == "constant N") {
      const auto Gamma = (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0))/
                         static_cast<T>(EquationData::Cp_Cv);

      return (parameters.p_bar/data.p_ref)*std::pow(BackgroundState::exner_constant_N<T>(p[dim - 1], data.L_ref,
                                                                                         parameters.N, parameters.T_bar),
                                                    static_cast<T>(1.0)/Gamma);
    }

    // Inverse of the scale height of the isothermal atmosphere
    const auto delta = static_cast<T>(EquationData::g)/(static_cast<T>(EquationData::R)*parameters.T_bar);

    return (parameters.p_bar/data.p_ref)*std::exp(-delta*p[dim - 1]*data.L_ref);
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
    Density(RunTimeParameters::IGWParameters& parameters_,
            RunTimeParameters::Data_Storage& data_);

    /**
     * Evaluation of the density
     * @param p point coordinate
     * @param component spatial component to be evaluated (unused)
     */
    virtual T value(const Point<dim, T>& p,
                    const unsigned       component = 0) const override;

  private:
    RunTimeParameters::IGWParameters& parameters; /*!< Parameters specific of the test case */
    RunTimeParameters::Data_Storage& data;        /*!< 'Global' parameters */
  };

  // Constructor which again relies on the 'Function' constructor.
  //
  template<unsigned dim, typename T>
  Density<dim, T>::Density(RunTimeParameters::IGWParameters& parameters_,
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

    const auto z = p[dim - 1]*data.L_ref;
    const auto r = (p[0]*data.L_ref - parameters.xc)/parameters.ac;

    // Background density (at z = 0)
    const auto rho_bar = parameters.p_bar/(static_cast<T>(EquationData::R)*parameters.T_bar);

    // Background state with constant buoyancy frequency and perturbation of the potential temperature
    if(parameters.background == "constant N") {
      const auto pi_bar      = BackgroundState::exner_constant_N<T>(p[dim - 1], data.L_ref, parameters.N, parameters.T_bar);
      const auto theta_bar   = BackgroundState::potential_temperature_constant_N<T>(p[dim - 1], data.L_ref,
                                                                                    parameters.N, parameters.T_bar);
      const auto theta_prime = parameters.dT*std::sin(numbers::PI*z/parameters.H)/(static_cast<T>(1.0) + r*r);

      return (rho_bar/data.rho_ref)*
             parameters.T_bar/(theta_bar + theta_prime)*std::pow(pi_bar, static_cast<T>(1.0)/
                                                                 (static_cast<T>(EquationData::Cp_Cv) - static_cast<T>(1.0)));
    }

    // Inverse of the scale height of the isothermal atmosphere
    const auto delta = static_cast<T>(EquationData::g)/(static_cast<T>(EquationData::R)*parameters.T_bar);

    const auto rho_0 = rho_bar*std::exp(-delta*z);

    // Perturbation of the temperature and corresponding perturbation of the density
    const auto T_prime = parameters.dT*std::exp(-r*r)*std::sin(numbers::PI*z/parameters.H);

    const auto rho_prime = -std::exp(-static_cast<T>(0.5)*delta*z)*rho_bar*T_prime/parameters.T_bar;

    return (rho_0 + rho_prime)/data.rho_ref;
  }

} // namespace IC_IGW

/**
 * @brief Inertia-gravity waves test case
 *
 * This type owns all inertia-gravity waves-specific configuration and hides the concrete
 * implementations from the solver.
 */
template<unsigned dim, typename T = double>
using IGWTestCase = TestCase<dim, T,
                             RunTimeParameters::IGWParameters,
                             IC_IGW::Density<dim, T>,
                             IC_IGW::Velocity<dim, T>,
                             IC_IGW::Pressure<dim, T>>;
