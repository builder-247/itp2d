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

#ifndef _TEST_RNG_HPP_
#define _TEST_RNG_HPP_

#include <utility>
#include <H5Cpp.h>
#include "tests_common.hpp"
#include "rng.hpp"

std::pair<double,double> mean_and_variance(std::vector<double> const& vec);

void write_sample(std::vector<double> const& vec, std::string filename);

#endif // _TEST_RNG_HPP_
