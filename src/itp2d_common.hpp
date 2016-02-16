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
 * Common type definitions, constants, and helper functions.
 */

#ifndef _ITP2D_COMMON_HPP_
#define _ITP2D_COMMON_HPP_

#include <limits>
#include <typeinfo>
#include <complex>
#include <fftw3.h>

#ifndef ITP2D_VERSION
#define ITP2D_VERSION "unknown"
#endif

extern const char version_string[];

__attribute__((unused)) const double pi = M_PI;
__attribute__((unused)) const double inf = std::numeric_limits<double>::infinity();
__attribute__((unused)) const double NaN = std::numeric_limits<double>::quiet_NaN();
__attribute__((unused)) const double machine_epsilon = std::numeric_limits<double>::epsilon();

/* Alias comp to the standard complex type.
 * NOTE: If your compiler stores std::complex<double> as anything other than
 * double[2], you are in big trouble. Luckily all compilers seem to follow this
 * de-facto standard.
 */
typedef std::complex<double> comp;

// Available boundary conditions
enum BoundaryType { Periodic, Dirichlet };

// Available orthonormalization algorithms
enum OrthoAlgorithm { Default, HighMem };

// Default FFTW flags
const unsigned int default_fftw_flags = FFTW_PATIENT;

// fftw_execute_dft for comp* arrays (why oh why can't FFTW use std::complex<double>?)
inline void fftw_execute_dft(const fftw_plan p, comp* in, comp* out) {
	fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}

// Shorthand for fftw_malloc
inline comp* malloc_comp(size_t N) {
	return reinterpret_cast<comp*>(fftw_malloc(N*sizeof(comp)));
}

// Simple rounding function
inline int round_to_int(double r) {
  return (int)((r > 0)? floor(r+0.5) : ceil(r-0.5));
}

#endif // _ITP2D_COMMON_HPP_
