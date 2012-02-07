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

/*
 * Unit tests for the StateSet class.
 */

#include "test_stateset.hpp"

// Check that orthonormalization works.
TEST(stateset, orthonormalization) {
	RNG rng(RNG::produce_random_seed());
	const DataLayout dl(16, 16, 1.0);
	StateSet states(8, dl, Default);
	states.init_to_gaussian_noise(rng);
	states.orthonormalize();
	ASSERT_LT(states.how_orthonormal(), 10*machine_epsilon);
}

TEST(stateset, orthonormalization_highmem) {
	RNG rng(RNG::produce_random_seed());
	const DataLayout dl(16, 16, 1.0);
	StateSet states(8, dl, HighMem);
	states.init_to_gaussian_noise(rng);
	states.orthonormalize();
	ASSERT_LT(states.how_orthonormal(), 10*machine_epsilon);
}
