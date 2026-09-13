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
 * Author: Giuseppe Orlando, 2026
 */
#pragma once

// Start by including the necessary header file
//
#include "../io/runtime_parameters.h"

// Forward declaration (so as not to include the header, not needed here
// and leave the possibility of other structures)
namespace RunTimeParameters {
  struct MountainData;
}

/**
 * @brief Auxiliary struct for the initial condition to implement
          the concrete instances of density, velocity, and pressure
 */
template<unsigned dim, typename T = double>
struct InitialCondition {
  std::unique_ptr<dealii::Function<dim, T>> density;
  std::unique_ptr<dealii::Function<dim, T>> velocity;
  std::unique_ptr<dealii::Function<dim, T>> pressure;
};

/**
 * @brief Abstract interface between the solver and a test case
 */
template<unsigned dim, typename T = double>
struct TestCaseBase {
  /**
   * Virtual destructor (class only for polymorphism)
   */
  virtual ~TestCaseBase() = default;

  InitialCondition<dim, T> ic; /*!< Auxiliary structure for concrete intial condition functions */

  /*!< Data for the (terrain-following) mesh mapping and the
       Rayleigh damping profiles. Defaults to "flat, undamped": a test
       case whose Parameters class does not opt in (see
       internal::has_mountain_data below) simply leaves this untouched. */
  std::unique_ptr<RunTimeParameters::MountainData> mountain_data;
};

namespace internal {
  // Detects whether 'Parameters' opts into supplying mountain/damping data,
  // i.e. exposes a mountain_data member, i.e. MountainData.
  // Same detection-idiom spirit as the is_constructible_v checks below,
  // just applied to a member function instead of a constructor. This is
  // how a Parameters class (e.g. MountainParameters) can hand 
  // their inputs without TestCase - or the solver - needing to
  // know it's "mountain" data.
  template<typename Parameters, typename = void>
  struct has_mountain_data: std::false_type {};
 
  template<typename Parameters>
  struct has_mountain_data<Parameters,
                           std::void_t<decltype(std::declval<Parameters&>().mountain_data)>>: std::true_type {};

  // Build a Function<dim,T>-derived 'Function' out of whichever subset of
  // (parameters, data) its constructor actually accepts: both, only one,
  // or neither. This is what allows density/velocity/pressure functions of
  // a test case to have different constructor needs without TestCase
  // having to know about it.
  //
  template<typename Parameters,
           class Function>
  std::unique_ptr<Function>
  make_ic_function(Parameters& parameters,
                   RunTimeParameters::Data_Storage& data) {
    if constexpr(std::is_constructible_v<Function, Parameters&, RunTimeParameters::Data_Storage&>) {
      return std::make_unique<Function>(parameters, data);
    }
    else if constexpr(std::is_constructible_v<Function, Parameters&>) {
      return std::make_unique<Function>(parameters);
    }
    else if constexpr(std::is_constructible_v<Function, RunTimeParameters::Data_Storage&>) {
      return std::make_unique<Function>(data);
    }
    else {
      return std::make_unique<Function>();
    }
  }
} // namespace internal

/**
 * @brief Generic test case, parametrized by its parameters type and by the
 *        concrete density/velocity/pressure function templates.
 *
 * Every concrete test case observed so far shares the same constructor
 * (parse a dedicated .prm file into a 'Parameters' object) and the same
 * setup() (build density/velocity/pressure out of that object). Rather
 * than repeating this boilerplate in each test case, a new test case only
 * has to provide:
 *   - a Parameters class exposing declare_parameters()/parse_parameters(),
 *   - Density/Velocity/Pressure function templates taking
 *     (Parameters&, RunTimeParameters::Data_Storage&) in their constructor,
 * and then declare a type alias, e.g.
 *
 *   template<unsigned dim, typename T = double>
 *   using NonHydrostatic3DTestCase = TestCase<dim, T,
                                               RunTimeParameters::MountainParameters,
                                               ICBC::Density<dim, T>,
                                               ICBC::Velocity<dim, T>,
                                               ICBC::Pressure<dim, T>>;
 *
 * A test case with genuinely different setup logic (e.g. needing extra
 * initialization not covered here) can still bypass this template and
 * derive from TestCaseBase directly.
 */
template<unsigned dim, typename T,
         class Parameters,
         class DensityFunction,
         class VelocityFunction,
         class PressureFunction>
class TestCase: public TestCaseBase<dim, T> {
public:
  /**
   * Class constructor: declares and parses the test-case-specific
   * parameter file into 'parameters'.
   * @param tc_param_file name of the parameter file
   * @param data 'global' parameters in case needed
   */
   explicit TestCase(const std::string& tc_param_file,
                     RunTimeParameters::Data_Storage& data);

protected:
  Parameters parameters; /*!< Auxiliary structure for test-case-dependent parameters */
};

// Class constructor and parsing of parameters (identical for every test case
// built on top of this template)
//
template<unsigned dim, typename T,
         class Parameters,
         class DensityFunction,
         class VelocityFunction,
         class PressureFunction>
TestCase<dim, T,
         Parameters,
         DensityFunction,
         VelocityFunction,
         PressureFunction>::
TestCase(const std::string& tc_param_file,
         RunTimeParameters::Data_Storage& data) {
  dealii::ParameterHandler prm;
  parameters.declare_parameters(prm);

  try {
    std::ifstream tc_file(tc_param_file);
    AssertThrow(tc_file, dealii::ExcFileNotOpen(tc_param_file));

    prm.parse_input(tc_file);
  }
  catch(const dealii::ExcFileNotOpen& exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Using default values for the parameters of the test case!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
  }

  parameters.parse_parameters(prm);

  if constexpr(internal::has_mountain_data<Parameters>::value) {
    this->mountain_data = std::make_unique<RunTimeParameters::MountainData>(parameters.mountain_data);
  }

  this->ic.density  = internal::make_ic_function<Parameters, DensityFunction>(parameters, data);
  this->ic.velocity = internal::make_ic_function<Parameters, VelocityFunction>(parameters, data);
  this->ic.pressure = internal::make_ic_function<Parameters, PressureFunction>(parameters, data);
}

// Printer of helper of parameters (identical for every test case)
//
template<class Parameters>
void print_help_tc() {
  dealii::ParameterHandler prm;
  Parameters parameters;

  parameters.declare_parameters(prm);

  prm.print_parameters(std::cout, dealii::ParameterHandler::Description);
}
