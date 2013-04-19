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
 * A test program for the EigenSolver class.
 */

#include "test_eigensolver.hpp"

TEST(eigensolver, eigenvalues_2x2) {
	// Test against some known eigenvalues of a 2x2 matrix.
	EigenSolver E(2);
	comp* matrix = new comp[4];
	matrix[0] = 1;
	matrix[1] = 2;
	matrix[2] = 2;
	matrix[3] = 4;
	E.solve(matrix);
	EXPECT_EQ(E.eigenvalue(0), 0);
	EXPECT_EQ(E.eigenvalue(1), 5);
	EXPECT_EQ(E.eigenvector(matrix, 0, 0), comp(-2/sqrt(5)));
	EXPECT_EQ(E.eigenvector(matrix, 0, 1), comp(1/sqrt(5)));
	EXPECT_EQ(E.eigenvector(matrix, 1, 0), comp(1/sqrt(5)));
	EXPECT_EQ(E.eigenvector(matrix, 1, 1), comp(2/sqrt(5)));
	delete[] matrix;
}

TEST(eigensolver, spectral_decomposition) {
	// Test that spectral decomposition works with the EigenSolver.
	const int N = 8;
	EigenSolver E(N);
	comp* matrix = new comp[N*N];
	comp* orig_matrix = new comp[N*N];
	comp* diff_matrix = new comp[N*N];
	for (int i=0; i<N; i++) {
		for (int j=0; j<i; j++)
			matrix[N*j+i] = comp(0, -abs(i-j));
		matrix[N*i+i] = comp(1, 0);
		for (int j=i+1; j<N; j++)
			matrix[N*j+i] = comp(0, abs(i-j));
	}
	memcpy(orig_matrix, matrix, N*N*sizeof(comp));
	E.solve(matrix);
	comp z;
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			z = 0;
			for (int n=0; n<N; n++) {
				z += E.eigenvector(matrix, n, i)*E.eigenvalue(n)*conj(E.eigenvector(matrix, n, j));
			}
			diff_matrix[N*j+i] = orig_matrix[N*j+i] - z;
		}
	}
	const double diff = cblas_dznrm2(N*N, reinterpret_cast<const double*>(diff_matrix), 1);
	delete[] matrix;
	delete[] orig_matrix;
	delete[] diff_matrix;
	EXPECT_LT(diff,  300*machine_epsilon);
}
