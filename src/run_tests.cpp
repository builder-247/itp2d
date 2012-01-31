/* Copyright 2012 Perttu Luukko

 * This file is part of itp2d.

 * itp2d is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.

 * itp2d is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.

 * You should have received a copy of the GNU General Public License along with
 * itp2d.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gtest/gtest.h>

// Include individual unit test suites

#include "test_rng.hpp"
#include "test_eigensolver.hpp"
#include "test_datalayout.hpp"
#include "test_statearray.hpp"
#include "test_state.hpp"
#include "test_transformer.hpp"
#include "test_stateset.hpp"
#include "test_operators.hpp"
#include "test_commandlineparser.hpp"
#include "test_potentialparser.hpp"
#include "test_laplacian.hpp"
#include "test_itp.hpp"
#include "test_hamiltonian.hpp"

// Global variable that determines whether unit tests dump internal data for deeper analysis
bool dump_data = false;

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  if (argc == 2 and std::string(argv[1]) == "--dump")
	  dump_data = true;
  return RUN_ALL_TESTS();
}
