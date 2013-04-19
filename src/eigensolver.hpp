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
 * A simple wrapper for LAPACK's ZHEEV solver
 */

#ifndef _EIGENSOLVER_HPP_
#define _EIGENSOLVER_HPP_

#include <cassert>
#include "itp2d_common.hpp"
#include "exceptions.hpp"

// MKL likes to use its own datatypes so we have to please it.

#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapack.h>
inline void my_zheev(char* jobz, char* uplo, const int* n, comp* a, const int*
		lda, double* w, comp* work, const int* lwork, double* rwork, int*
		info) {
	zheev(jobz, uplo, const_cast<int*>(n), reinterpret_cast<MKL_Complex16*>(a), const_cast<int*>(lda), w,
			reinterpret_cast<MKL_Complex16*>(work), const_cast<int*>(lwork), rwork, info);
}
#else
extern "C" {
#include <cblas.h>
extern void zheev_(char* jobz, char* uplo, const int* n, comp* a, const int*
		lda, double* w, comp* work, const int* lwork, double* rwork, int*
		info);
}
inline void my_zheev(char* jobz, char* uplo, const int* n, comp* a, const int*
		lda, double* w, comp* work, const int* lwork, double* rwork, int*
		info) {
	zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
#endif

// A solver for eigenvalues and -vectors of NxN complex Hermitian matrices

class EigenSolver {
	public:
		EigenSolver(size_t N);
		~EigenSolver();
		inline comp const& eigenvector(comp const* input_matrix, size_t n, size_t i) const { return input_matrix[n*size+i]; } // i:th element of n:th eigenvector
		inline void scale_eigenvector(comp* input_matrix, size_t n, double value) const;
		inline double const& eigenvalue(size_t n) const { return evals[n]; }
		inline void solve(comp* input_matrix); // Note: Input data must be specified in column-major (FORTRAN) order! Also note that this destroys the matrix.
	private:
		const int size;
		int lwork_size;
		int info;
		double* evals;
		comp* lwork;
		double* rwork;
		static char DoIWantVectors[2];
		static char UpperOrLower[2];
};

inline void EigenSolver::scale_eigenvector(comp* input_matrix, size_t n, double value) const {
	cblas_zdscal(size, value, reinterpret_cast<double*>(input_matrix+n*size), 1);
}

inline void EigenSolver::solve(comp* input_matrix) {
	my_zheev(DoIWantVectors, UpperOrLower, &size, input_matrix, &size, evals, lwork, &lwork_size, rwork, &info);
	if (info != 0)
		throw EigensolverError(info);
}

#endif // _EIGENSOLVER_HPP_
