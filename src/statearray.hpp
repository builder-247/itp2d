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
 * A class representing an array of State objects on a common DataLayout. This
 * is the underlying low-level object under StateSet. Also used to hold
 * temporary workspace data.
 */

#ifndef _STATEARRAY_HPP_
#define _STATEARRAY_HPP_

#include "state.hpp"

class StateArray {
	public:
		StateArray(size_t N, DataLayout const& dl); // Array of N states with a common DataLayout
		StateArray(size_t N, DataLayout const& dl, comp* ptr);	// Use the provided pointer for storage
		// Constructors for creating slices of existing StateArray objects.
		StateArray(StateArray& arr, size_t start_index);	// A slice of all the states starting from a certain index
		StateArray(StateArray& arr, size_t start_index, size_t len); // ...or only len states starting from a certain index
		~StateArray();
		inline comp* get_dataptr() { return memptr; }
		inline State& operator[](size_t n) { assert(n<N); return *(ptrarray[n]); }
		inline State const& operator[](size_t n) const { assert(n<N); return *(ptrarray[n]); }
		inline size_t size() const { return N; }
		DataLayout const& datalayout;
	private:
		static State** allocate_ptrarray(size_t N, DataLayout const& dl, comp* const memptr);
		const size_t N;
		const bool is_slice;
		const bool handle_own_memory;
		comp* const memptr;
		State* const* const ptrarray;
};

#endif // _STATEARRAY_HPP_
