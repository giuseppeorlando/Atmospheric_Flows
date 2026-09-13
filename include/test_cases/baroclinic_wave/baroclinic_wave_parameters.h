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

#include <deal.II/base/parameter_handler.h>

// @sect{Run-time parameters}
//
// Since our method has several parameters that can be fine-tuned we put them
// into an external file, so that they can be determined at run-time.
//
namespace RunTimeParameters {
  using namespace dealii;
  
  /**
   * @brief Parser for the parameters related to mountain test cases
   *
   * These do not belong in the global parser Data_Storage.
   */
  class BaroclinicWaveParameters {
  public:
    /**
     * Class constructor
     */
    BaroclinicWaveParameters() = default;

    /**
     * Declare the parameters for mountain test cases
     */
    void declare_parameters(ParameterHandler& prm);

    /**
     * Parse the parameters for mountain test cases
     */
    void parse_parameters(ParameterHandler& prm);

    /*--- Now declare the parameters to be parsed ---*/
    // Physical parameters
    double u_bar;   /*!< Reference background velocity */
    double p_bar;   /*!< Reference background pressure */
    double T_bar;   /*!< Reference background temperature */

    double lapse_rate;
    double b;
    double Omega0;
    double phi0;
    double up;
    double Lp;
    double xc;
    double yc;
  };

  // Declare the parameters for mountain test cases
  //
  void BaroclinicWaveParameters::declare_parameters(ParameterHandler& prm) {
    // Physical parameters
    prm.enter_subsection("Physical data");
    {
      prm.declare_entry("u_bar",
                        "1.0",
                        Patterns::Double(),
                        "The (horizontal) background.");
      prm.declare_entry("p_bar",
                        "1.0",
                        Patterns::Double(0.0),
                        "The background pressure (at z = 0).");
      prm.declare_entry("T_bar",
                        "1.0",
                        Patterns::Double(0.0),
                        "The background temperature (at z = 0).");
      
      prm.declare_entry("lapse_rate",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("b",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("Omega0",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("phi0",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("up",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("Lp",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("xc",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");
      prm.declare_entry("yc",
                        "1.0",
                        Patterns::Double(0.0),
                        "boh.");                  
    }
    prm.leave_subsection();
  }

  // Parse the parameters for mountain test cases
  //
  void BaroclinicWaveParameters::parse_parameters(ParameterHandler& prm) {
    // Physical parameters
    prm.enter_subsection("Physical data");
    {
      u_bar = prm.get_double("u_bar");
      p_bar = prm.get_double("p_bar");
      T_bar = prm.get_double("T_bar");

      lapse_rate = prm.get_double("lapse_rate");
      b          = prm.get_double("b");
      Omega0     = prm.get_double("Omega0");
      phi0       = prm.get_double("phi0");
      up         = prm.get_double("up");
      Lp         = prm.get_double("Lp");
      xc         = prm.get_double("xc");
      yc         = prm.get_double("yc");
    }
    prm.leave_subsection();
  }

} // namespace RunTimeParameters
