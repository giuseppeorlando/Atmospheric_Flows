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

#include <deal.II/base/parameter_handler.h>

#include <string>

// @sect{Run-time parameters}
//
// Since our method has several parameters that can be fine-tuned we put them
// into an external file, so that they can be determined at run-time.
//
namespace RunTimeParameters {
  using namespace dealii;

  /**
   * @brief Parser for the parameters related to the inertia-gravity waves test case
   *
   * These do not belong in the global parser Data_Storage.
   */
  class IGWParameters {
  public:
    /**
     * Class constructor
     */
    IGWParameters() = default;

    /**
     * Declare the parameters for the inertia-gravity waves test case
     */
    void declare_parameters(ParameterHandler& prm);

    /**
     * Parse the parameters for the inertia-gravity waves test case
     */
    void parse_parameters(ParameterHandler& prm);

    /*--- Now declare the parameters to be parsed ---*/
    // Parameters of the perturbation
    double dT; /*!< Amplitude of the perturbation of the (potential) temperature (0 for a background state at rest) */
    double H;  /*!< Vertical scale of the perturbation (typically the height of the domain) */
    double xc; /*!< x-center of the perturbation */
    double ac; /*!< Semi-width of the perturbation */

    // Background state
    std::string background; /*!< Type of background state ('isothermal' or 'constant N') */
    double N;               /*!< Buoyancy frequency (only for 'constant N', N = 0 for a neutral atmosphere) */
    double u_bar;           /*!< Reference background velocity */
    double p_bar;           /*!< Reference background pressure */
    double T_bar;           /*!< Reference background temperature ('isothermal')
                                 or potential temperature at z = 0 ('constant N') */
  };

  // Declare the parameters for the inertia-gravity waves test case
  //
  void IGWParameters::declare_parameters(ParameterHandler& prm) {
    prm.enter_subsection("Physical data");
    {
      prm.declare_entry("dT",
                        "0.01",
                        Patterns::Double(),
                        "The amplitude of the perturbation of the temperature ('isothermal') or of the potential "
                        "temperature ('constant N'). Set it to 0 for a background state at rest.");
      prm.declare_entry("H",
                        "10000.0",
                        Patterns::Double(0.0),
                        "The vertical scale of the perturbation (typically the height of the domain).");
      prm.declare_entry("xc",
                        "0.0",
                        Patterns::Double(),
                        "The x-center of the perturbation.");
      prm.declare_entry("ac",
                        "1.0",
                        Patterns::Double(0.0),
                        "The semi-width of the perturbation.");

      prm.declare_entry("background",
                        "isothermal",
                        Patterns::Selection("isothermal|constant N"),
                        "The type of background state: 'isothermal' or 'constant N' (constant buoyancy frequency).");
      prm.declare_entry("N",
                        "0.01",
                        Patterns::Double(0.0),
                        "The buoyancy frequency (only for 'constant N', N = 0 for a neutral atmosphere).");
      prm.declare_entry("u_bar",
                        "0.0",
                        Patterns::Double(),
                        "The (horizontal) background velocity.");
      prm.declare_entry("p_bar",
                        "100000.0",
                        Patterns::Double(0.0),
                        "The background pressure (at z = 0).");
      prm.declare_entry("T_bar",
                        "250.0",
                        Patterns::Double(0.0),
                        "The background temperature ('isothermal') or the potential temperature at z = 0 ('constant N').");
    }
    prm.leave_subsection();
  }

  // Parse the parameters for the inertia-gravity waves test case
  //
  void IGWParameters::parse_parameters(ParameterHandler& prm) {
    prm.enter_subsection("Physical data");
    {
      dT = prm.get_double("dT");
      H  = prm.get_double("H");
      xc = prm.get_double("xc");
      ac = prm.get_double("ac");

      background = prm.get("background");
      N          = prm.get_double("N");
      u_bar      = prm.get_double("u_bar");
      p_bar = prm.get_double("p_bar");
      T_bar = prm.get_double("T_bar");
    }
    prm.leave_subsection();
  }

} // namespace RunTimeParameters
