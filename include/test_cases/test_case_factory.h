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

// Start by including the necessary header files
//
#include "mountain/ic_3D_nonhydrostatic.h"
#include "mountain/Rayleigh_damping.h"

#include <filesystem>
namespace fs = std::filesystem;

/**
 * @brief Construct a concrete test case from its user-facing name
 *
 * Keep this as a simple if-chain initially. It can later become a registry
 * without changing the solver-side interface.
 */
template<unsigned dim, typename T = double, typename... Args>
std::unique_ptr<TestCaseBase<dim, T>>
make_test_case(const std::string& tc_name,
               const std::string& tc_param_file,
               Args &&...args) {
  if(tc_name == "3D non-hydrostatic mountain") {
    return std::make_unique<NonHydrostatic3DTestCase<dim, T>>(tc_param_file, std::forward<Args>(args)...);
  }

  AssertThrow(false,
              dealii::ExcMessage("Unknown test case: " + tc_name));
}

/**
 * @brief Print helper of a concrete test case from its user-facing name
 *
 * Keep this as a simple if-chain initially. It can later become a registry
 * without changing the solver-side interface.
 */
void print_help_test_case(const std::string& tc_name) {
  if(tc_name == "3D non-hydrostatic mountain") {
    return print_help_tc<RunTimeParameters::MountainParameters>();
  }

  AssertThrow(false,
              dealii::ExcMessage("Unknown test case: " + tc_name));
}
