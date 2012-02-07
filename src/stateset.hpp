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
 * A class representing a set of State objects. This is the object that is worked on by the ITP algorithm.
 * Orthonormalization of states is also defined here.
 */

#ifndef _STATESET_HPP_
#define _STATESET_HPP_

#include <vector>
#include <utility>
#include <map>
#include <exception>
#include <cmath>
#include <omp.h>
#include "H5Cpp.h"

#ifdef USE_MKL
#include <mkl_cblas.h>
#else
extern "C" {
#include <cblas.h>
}
#endif

#include "itp2d_common.hpp"
#include "exceptions.hpp"
#include "state.hpp"
#include "statearray.hpp"
#include "timer.hpp"
#include "rng.hpp"
#include "eigensolver.hpp"
#include "parameters.hpp"

class StateSet {
	public:
		StateSet(size_t N, DataLayout const& dl, OrthoAlgorithm algo = Default);
		~StateSet();
		// Initializing
		void init(Parameters const& params, RNG& rng);
		void init(comp (*initfunc)(size_t n, double x, double y)); // Initialize from function
		void init_from_datafile(std::string filename);
		void init_to_gaussian_noise(RNG& rng);
		// Simple getters & setters
		inline State& operator[](size_t n) const { return (*state_array)[n]; }
		inline size_t get_num_states() const { return N; }
		inline void set_timestep_converged(size_t n, bool val=true);
		inline bool is_timestep_converged(size_t n) const { return timestep_converged[n]; }
		inline size_t get_num_timestep_converged() const { return how_many_timestep_converged; }
		inline void set_finally_converged(size_t n, bool val=true);
		inline bool is_finally_converged(size_t n) const { return finally_converged[n]; }
		inline size_t get_num_finally_converged() const { return how_many_finally_converged; }
		// Arithmetic
		inline comp dot(size_t i, size_t j) const;
		// Orthonormalizing
		void orthonormalize() throw(std::exception);
		bool is_orthonormal(double epsilon = 1e-5) const;
		double how_orthonormal() const;
		// Timing
		double get_ortho_time() const { return ortho_time; };
		double get_dot_time() const { return dot_time; };
		double get_eigensolve_time() const { return eigensolve_time; };
		double get_lincomb_time() const { return lincomb_time; };
		// Public reference to the underlying data layout
		DataLayout const& datalayout;
	private:
		const size_t N;
		const OrthoAlgorithm ortho_algorithm;
		// Normally all operations are done as much in-place as possible to
		// conserve memory. However, with the HighMem OrthoAlgorithm operations
		// are done out-of-place, and for that we need to store the data
		// *twice* (once for the input of an operation, and once for the
		// output). This allows for possibly better performance in the expense
		// of roughly doubled memory requirement.
		StateArray* state_array;
		StateArray* other_state_array;
		comp* dataptr1;
		comp* dataptr2;
		StateArray* statearrayptr1;
		StateArray* statearrayptr2;
		EigenSolver ESolver;
		comp* overlapmatrix;
		std::vector<comp> tempstate;
		std::vector<bool> timestep_converged;
		std::vector<bool> finally_converged;
		size_t how_many_timestep_converged;
		size_t how_many_finally_converged;
		inline comp& data(size_t n, size_t x, size_t y) { return (*state_array)[n](x,y); }
		inline void switch_state_arrays();
		// For timing
		double ortho_time, dot_time, eigensolve_time, lincomb_time;
		Timer ortho_timer, dot_timer, eigensolve_timer, lincomb_timer;
};

// Simple getters & setters

inline void StateSet::set_finally_converged(size_t n, bool val) {
	if (finally_converged[n] == val)
		return;
	if (val) {
		finally_converged[n] = true;
		how_many_finally_converged++;
	}
	else {
		finally_converged[n] = false;
		how_many_finally_converged--;
	}
}

inline void StateSet::set_timestep_converged(size_t n, bool val) {
	if (timestep_converged[n] == val)
		return;
	if (val) {
		timestep_converged[n] = true;
		how_many_timestep_converged++;
	}
	else {
		timestep_converged[n] = false;
		how_many_timestep_converged--;
	}
}

// Arithmetic

inline comp StateSet::dot(size_t i, size_t j) const {
	return (*state_array)[i].dot((*state_array)[j]);
}

inline void StateSet::switch_state_arrays() {
	other_state_array = state_array;
	state_array = (state_array == statearrayptr2)? statearrayptr1 : statearrayptr2;
}

#endif // _STATESET_HPP_
