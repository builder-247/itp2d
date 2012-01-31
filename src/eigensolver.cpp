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
 * A simple wrapper for LAPACK's ZHEEV eigenvalue solver
 */

#include "eigensolver.hpp"

char EigenSolver::DoIWantVectors[] = "V"; // Yes I want eigenvectors
char EigenSolver::UpperOrLower[] = "U"; // Use upper triangular part

EigenSolver::EigenSolver(size_t N) : size(static_cast<int>(N)) {
	// Allocate temporary space required by ZHEEV
	// These space requirements of these two arrays are always known
	evals = new double[N];
	rwork = new double[3*N-2];
	// For lwork we need to first query the optimal size by calling ZHEEV once
	// with a NULL argument. This is stupid but this is how LAPACK does it.
	comp temp = 0;
	lwork_size = -1;
	lwork = &temp;
	solve(NULL);	
	assert(info == 0);	// The NULL-run should never fail but let's be doubly sure.
	// Now work[0] should contain the optimal size.
	lwork_size = static_cast<int>(std::real(lwork[0]));
	lwork = new comp[lwork_size];
}

EigenSolver::~EigenSolver() {
	delete[] evals;
	delete[] lwork;
	delete[] rwork;
}
