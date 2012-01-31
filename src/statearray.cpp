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

#include "statearray.hpp"

StateArray::StateArray(size_t arg_N, DataLayout const& dl) :
		datalayout(dl),
		N(arg_N),
		is_slice(false),
		handle_own_memory(true),
		memptr(malloc_comp(N*datalayout.N)),
		ptrarray(allocate_ptrarray(N, datalayout, memptr)) {}

StateArray::StateArray(size_t arg_N, DataLayout const& dl, comp* ptr) :
		datalayout(dl),
		N(arg_N),
		is_slice(false),
		handle_own_memory(false),
		memptr(ptr),
		ptrarray(allocate_ptrarray(N, datalayout, memptr)) {}

StateArray::StateArray(StateArray& arr, size_t start_index) :
		datalayout(arr.datalayout),
		N(arr.N-start_index),
		is_slice(true),
		handle_own_memory(false),
		memptr((N==0) ? NULL : arr.memptr + start_index*datalayout.N),
		ptrarray(arr.ptrarray + start_index) {
	assert(start_index <= arr.N);
}

StateArray::StateArray(StateArray& arr, size_t start_index, size_t len) :
		datalayout(arr.datalayout),
		N(len),
		is_slice(true),
		handle_own_memory(false),
		memptr((N==0) ? NULL : arr.memptr + start_index*datalayout.N),
		ptrarray(arr.ptrarray + start_index) {
	assert(start_index+len <= arr.N);
}

StateArray::~StateArray() {
	if (not is_slice) {
		for (size_t i=0; i<N; i++) {
			delete ptrarray[i];
		}
		delete[] ptrarray;
		if (handle_own_memory) {
			fftw_free(memptr);
		}
	}
}

State** StateArray::allocate_ptrarray(size_t N, DataLayout const& dl, comp* const memptr) { 
	State** ptrarray = new State*[N];
	for (size_t i=0; i<N; i++) {
		ptrarray[i] = new State(dl, memptr+i*dl.N);
	}
	return ptrarray;
}
