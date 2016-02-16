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
 * A simple class for storing the low-level layout of 2D data -- the x- and
 * y-lengths, the total number of grid points and the grid spacing (dx).
 */

#ifndef _DATALAYOUT_HPP_
#define _DATALAYOUT_HPP_

#include "itp2d_common.hpp"
#include "exceptions.hpp"

class DataLayout {
	public:
		DataLayout(size_t sx, size_t sy, double dx);
		~DataLayout();
		// Helper functions for getting values from a "flattened" 1D array, using this 2D data layout
		template<typename Type> inline Type& value(Type* array, size_t x, size_t y) const;
		template<typename Type> inline Type const& value(Type const* array, size_t x, size_t y) const;
		inline double const& get_posx(size_t x) const { return posx[x]; }
		inline double const& get_posy(size_t y) const { return posy[y]; }
		inline size_t get_x_index(double x) const { return nearest_index(x, sizex, dx); }
		inline size_t get_y_index(double y) const { return nearest_index(y, sizey, dx); }
		// const data members
		const size_t sizex;
		const size_t sizey;
		const size_t N;
		const double dx;
		const double lenx;
		const double leny;
	private:
		// these arrays store the positions of individual grid points
		double* const posx;	
		double* const posy;
		// private helper function for get_?_index
		static inline size_t nearest_index(double x, size_t s, double dx);
};

// Free functions for comparison testing
inline bool operator==(DataLayout const& lhs, DataLayout const& rhs) {
	if (&lhs == &rhs)
		return true;
	if (lhs.sizex == rhs.sizex and lhs.sizey == rhs.sizey and lhs.dx == rhs.dx)
		return true;
	return false;
}

inline bool operator!=(DataLayout const& lhs, DataLayout const& rhs) {
	return !(lhs == rhs);
}

template<typename Type>
inline Type& DataLayout::value(Type* array, size_t x, size_t y) const {
	return array[y*sizex+x];
}

template<typename Type>
inline Type const& DataLayout::value(Type const* array, size_t x, size_t y) const {
	return array[y*sizex+x];
}

inline size_t DataLayout::nearest_index(double x, size_t s, double d) {
	int i = round_to_int(x/d + (static_cast<double>(s)-1)/2);
	if (i < 0) {
		throw GeneralError("Overflow in DataLayout::nearest_index, returned index negative");
	}
	if (i >= static_cast<int>(s)) {
		throw GeneralError("Overflow in DataLayout::nearest_index, returned index too large");
	}
	return i;
}

#endif // _DATALAYOUT_HPP_
