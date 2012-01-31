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

#include "state.hpp"

// Free functions for comparison testing
bool operator==(State const& lhs, State const& rhs) {
	if (&lhs == &rhs)
		return true;
	if (lhs.datalayout != rhs.datalayout)
		return false;
	comp const* lhs_data = lhs.data_ptr();
	comp const* rhs_data = rhs.data_ptr();
	for (size_t n=0; n<lhs.datalayout.N; n++)
		if (lhs_data[n] != rhs_data[n])
			return false;
	return true;
}

bool operator!=(State const& lhs, State const& rhs) {
	return !(lhs == rhs);
}

// Distance functions for comparing two states based on either the root mean
// square distance or maximum distance.
double rms_distance(State const& lhs, State const& rhs) {
	if (lhs.datalayout != rhs.datalayout)
		return inf;
	comp const* lhs_data = lhs.data_ptr();
	comp const* rhs_data = rhs.data_ptr();
	double rmssum = 0;
	for (size_t n=0; n<lhs.datalayout.N; n++)
		rmssum += norm(lhs_data[n]-rhs_data[n]);
	rmssum /= static_cast<double>(lhs.datalayout.N);
	return sqrt(rmssum);
}

double max_distance(State const& lhs, State const& rhs) {
	if (lhs.datalayout != rhs.datalayout)
		return inf;
	comp const* lhs_data = lhs.data_ptr();
	comp const* rhs_data = rhs.data_ptr();
	double max = 0;
	for (size_t n=0; n<lhs.datalayout.N; n++) {
		const double dist = abs(lhs_data[n]-rhs_data[n]);
		if (dist > max)
			max = dist;
	}
	return max;
}

// Overloaded global operator for easy printing of State data

std::ostream& operator<<(std::ostream& stream, const State& state) {
	for (size_t y=0; y<state.datalayout.sizey; y++) {
		for (size_t x=0; x<state.datalayout.sizex-1; x++) {
			stream << state(x,y) << " ";
		}
		stream << state(state.datalayout.sizey-1, y);
		stream << std::endl;
	}
	return stream;
}

// Constructors & Destructors

State::State(DataLayout const& lay) :
		datalayout(lay),
		handle_own_memory(true),
		memptr(malloc_comp(datalayout.N)) {
}

State::State(const State& other) :
		datalayout(other.datalayout),
		handle_own_memory(true),
		memptr(malloc_comp(datalayout.N)) {
	assert(memptr != other.memptr);
	memcpy(memptr, other.memptr, datalayout.N*sizeof(comp));
}

State::State(DataLayout const& lay, comp (*initfunc)(double, double)) :
		datalayout(lay),
		handle_own_memory(true),
		memptr(malloc_comp(datalayout.N)) {
	set_by_func(initfunc);
}

State::State(DataLayout const& lay, comp* ptr) :
		datalayout(lay),
		handle_own_memory(false),
		memptr(ptr) {}

State::State(const State& other, comp* ptr) :
		datalayout(other.datalayout),
		handle_own_memory(false),
		memptr(ptr) {
	assert(memptr != other.memptr);
	memcpy(memptr, other.memptr, datalayout.N*sizeof(comp));
}

State::State(DataLayout const& lay, comp (*initfunc)(double, double), comp* ptr) :
		datalayout(lay),
		handle_own_memory(false),
		memptr(ptr) {
	set_by_func(initfunc);
}

State::~State() {
	if (handle_own_memory) {
		fftw_free(memptr);
	}
}
