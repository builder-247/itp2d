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

#include "rng.hpp"

RNG::RNG(unsigned long int s) : seed(s), uniform_rng(base_rng_type(seed), uniform_distribution_type()), normal_distribution() {}

RNG::RNG() : seed(produce_random_seed()), uniform_rng(base_rng_type(seed), uniform_distribution_type()), normal_distribution() {}

unsigned long int RNG::produce_random_seed() {
	timeval time_now;
	gettimeofday(&time_now, NULL);
	return static_cast<unsigned long int>(time_now.tv_usec);
}
