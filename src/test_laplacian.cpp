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
 * Unit tests for the Laplacian operator. This mainly tests whether the FFT
 * operations work correctly.
 */

#include "test_laplacian.hpp"

namespace test_laplacian_reference {
	comp initfunc(double x, double y) {
		return comp(exp(-(x*x + y*y)), 0);
	}

	comp referencefunc(double x, double y) {
		const double r2 = x*x + y*y;
		const double z = 4*(r2-1)*exp(-r2);
		return comp(z, 0);
	}
}

TEST(laplacian, laplacian_of_gaussian) {
	const size_t N = 42;
	const double len = 12;
	const double dx = len/N;
	const DataLayout dl(N, N, dx);
	const Transformer tr(dl);
	State G(dl, test_laplacian_reference::initfunc);
	StateArray work(0, dl);
	Laplacian Lapl(tr);
	Lapl(G, work);
	State R(dl, test_laplacian_reference::referencefunc);
	const double rmsdist = rms_distance(G,R);
	const double maxdist = max_distance(G,R);
	ASSERT_LT(rmsdist, 100*machine_epsilon);
	ASSERT_LT(maxdist, 1000*machine_epsilon);
	if (dump_data) {
		Datafile file("data/test_laplacian.h5", dl, true);
		file.write_state(0, 0, G);
		file.write_state(0, 1, R);
	}
}
