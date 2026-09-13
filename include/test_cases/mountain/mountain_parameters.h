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
   * @brief Whatever the (terrain-following) mesh mapping and the
   *        Rayleigh damping profiles need.
   *
   * This is deliberately configuration-agnostic: it knows nothing about
   * mountains. A test case that has no mountain and no sponge layers 
   * does not need to know this struct exists at all; 
   * it simply does not override it (see TestCaseBase::mountain_data in test_case.h).
   *
   * A test case that does need it (e.g. a mountain) exposes it by giving its
   * Parameters class a `MountainData get_mountain_data() const` method; see
   * RunTimeParameters::MountainParameters::get_mountain_data() for the
   * concrete example. TestCase<...> picks it up automatically if present
   * (again, see test_case.h).
   */
struct MountainData {
  // Physical parameters
  double h, xc, yc, ac;
 
  // Rayleigh damping parameters
  double z_start, lambda_z;
  double x_start_left, lambda_x_left;
  double x_start_right, lambda_x_right;
  double y_start_left, lambda_y_left;
  double y_start_right, lambda_y_right;
};

  /**
   * @brief Parser for the parameters related to mountain test cases
   *
   * These do not belong in the global parser Data_Storage.
   */
  class MountainParameters {
  public:
    /**
     * Class constructor
     */
    MountainParameters() = default;

    /**
     * Declare the parameters for mountain test cases
     */
    void declare_parameters(ParameterHandler& prm);

    /**
     * Parse the parameters for mountain test cases
     */
    void parse_parameters(ParameterHandler& prm);

    /*--- Now declare the parameters to be parsed ---*/
    MountainData mountain_data; /*!< Mountain and Rayleigh damping parameters */

    // Other physical parameters
    double N; /*!< Buoyancy frequency */

    double u_bar;   /*!< Reference background velocity */
    double p_bar;   /*!< Reference background pressure */
    double T_bar;   /*!< Reference background temperature */
    double rho_bar; /*!< Reference background density */
  };

  // Declare the parameters for mountain test cases
  //
  void MountainParameters::declare_parameters(ParameterHandler& prm) {
    // Physical parameters
    prm.enter_subsection("Physical data");
    {
      prm.declare_entry("h",
                        "1.0",
                        Patterns::Double(0.0),
                        "The mountain height.");
      prm.declare_entry("xc",
                        "1.0",
                        Patterns::Double(),
                        "The x-center of the mountain.");
      prm.declare_entry("yc",
                        "1.0",
                        Patterns::Double(),
                        "The y-center of the mountain.");
      prm.declare_entry("ac",
                        "1.0",
                        Patterns::Double(0.0),
                        "The semi-width of the mountain.");

      prm.declare_entry("N",
                        "0.01",
                        Patterns::Double(0.0),
                        "The buoyancy frequency.");

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
      prm.declare_entry("rho_bar",
                        "1.0",
                        Patterns::Double(0.0),
                        "The background density (at z = 0).");
    }
    prm.leave_subsection();

    // Rayleigh damping parameters
    prm.enter_subsection("Rayleigh damping data");
    {
      prm.declare_entry("z_start",
                        "1.0",
                        Patterns::Double(),
                        "Start of Rayleigh damping for top boundary.");
      prm.declare_entry("lambda_z",
                        "1.0",
                        Patterns::Double(0.0),
                        "Intensity of Rayleigh damping for top boundary.");
      prm.declare_entry("x_start_left",
                        "0.0",
                        Patterns::Double(),
                        "Start of Rayleigh damping for left boundary.");
      prm.declare_entry("lambda_x_left",
                        "1.0",
                        Patterns::Double(0.0),
                        "Intensity of Rayleigh damping for left boundary.");
      prm.declare_entry("x_start_right",
                        "1.0",
                        Patterns::Double(0.0),
                        "Start of Rayleigh damping for right boundary.");
      prm.declare_entry("lambda_x_right",
                        "1.0",
                        Patterns::Double(),
                        "Intensity of Rayleigh damping for right boundary.");
      prm.declare_entry("y_start_left",
                        "0.0",
                        Patterns::Double(0.0),
                        "Start of Rayleigh damping for y left boundary.");
      prm.declare_entry("lambda_y_left",
                        "1.0",
                        Patterns::Double(0.0),
                        "Intensity of Rayleigh damping for y left boundary.");
      prm.declare_entry("y_start_right",
                        "1.0",
                        Patterns::Double(0.0),
                        "Start of Rayleigh damping for y right boundary.");
      prm.declare_entry("lambda_y_right",
                        "1.0",
                        Patterns::Double(),
                        "Intensity of Rayleigh damping for y right boundary.");
    }
    prm.leave_subsection();
  }

  // Parse the parameters for mountain test cases
  //
  void MountainParameters::parse_parameters(ParameterHandler& prm) {
    // Physical parameters
    prm.enter_subsection("Physical data");
    {
      mountain_data.h  = prm.get_double("h");
      mountain_data.xc = prm.get_double("xc");
      mountain_data.yc = prm.get_double("yc");
      mountain_data.ac = prm.get_double("ac");

      N = prm.get_double("N");

      u_bar   = prm.get_double("u_bar");
      p_bar   = prm.get_double("p_bar");
      T_bar   = prm.get_double("T_bar");
      rho_bar = prm.get_double("rho_bar");
    }
    prm.leave_subsection();

    // Rayleigh damping parameters
    prm.enter_subsection("Rayleigh damping data");
    {
      mountain_data.z_start        = prm.get_double("z_start");
      mountain_data.lambda_z       = prm.get_double("lambda_z");
      mountain_data.x_start_left   = prm.get_double("x_start_left");
      mountain_data.lambda_x_left  = prm.get_double("lambda_x_left");
      mountain_data.x_start_right  = prm.get_double("x_start_right");
      mountain_data.lambda_x_right = prm.get_double("lambda_x_right");
      mountain_data.y_start_left   = prm.get_double("y_start_left");
      mountain_data.lambda_y_left  = prm.get_double("lambda_y_left");
      mountain_data.y_start_right  = prm.get_double("y_start_right");
      mountain_data.lambda_y_right = prm.get_double("lambda_y_right");
    }
    prm.leave_subsection();
  }

} // namespace RunTimeParameters
