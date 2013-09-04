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
 * A class representing a quantum wave function in a discrete 2D grid. This is
 * essentially a big array of complex numbers, with the gridded structure given
 * by the DataLayout class, and with basic arithmetic functions added for
 * ease-of-use.
 */

#ifndef _STATES_HPP_
#define _STATES_HPP_

#include <cassert>
#include <cstring>
#include "itp2d_common.hpp"
#include "datalayout.hpp"
#include "transformer.hpp"

#ifdef USE_MKL
#include <mkl_cblas.h>
#else
extern "C" {
#include <cblas.h>
}
#endif

class State {
	public:
		State(DataLayout const& lay);
		State(const State& other);
		State(DataLayout const& lay, comp (*initfunc)(double, double));
		State(DataLayout const& lay, comp* ptr);
		State(const State& other, comp* ptr);
		State(DataLayout const& lay, comp (*initfunc)(double, double), comp* ptr);
		~State();
		// Basic operations
		void operator=(State const& other);
		inline void zero();
		inline void normalize(double target_norm = 1.0);
		// Getters & Setters
		inline comp& operator()(size_t x, size_t y) { return datalayout.value(memptr, x, y); }
		inline comp const& operator()(size_t x, size_t y) const { return datalayout.value(memptr, x, y); }
		inline comp const* data_ptr() const { return memptr; }
		void set_by_func(comp (*initfunc)(double, double));
		// Arithmetic
		inline State& operator+=(const State& other);
		inline State& operator-=(const State& other);
		inline State& operator*=(const double& other);
		inline State& operator/=(const double& other);
		inline State& operator*=(const comp& other);
		inline State& operator/=(const comp& other);
		template<typename Type> inline void pointwise_multiply(Type const* values);
		template<typename Type> inline void pointwise_divide(Type const* values);
		template<typename Type> inline void pointwise_multiply_imaginary_shiftx(Type const* values);
		template<typename Type> inline void pointwise_multiply_y(Type const* values);
		template<typename Type> inline void pointwise_multiply_and_add(Type const* values, State const& addstate);
		inline comp dot(const State& other) const;
		inline double norm() const;
		// DFT, DST and DCT operations
		inline void transform(Transform trans, Transformer const& tr) {
			assert(tr.datalayout == datalayout);
			tr.transform(memptr, trans);
		}
		// Public reference to underlying DataLayout
		DataLayout const& datalayout;
	private:
		const bool handle_own_memory;
		comp* const memptr;
};

// Free functions for comparison testing
bool operator==(State const& lhs, State const& rhs);
bool operator!=(State const& lhs, State const& rhs);
double rms_distance(State const& lhs, State const& rhs);
double max_distance(State const& lhs, State const& rhs);

// Overloaded global operator for easy printing of State data

std::ostream& operator<<(std::ostream& stream, const State& state);

/*
 *		Arithmetic functions
 */

// Basic operations

inline void State::operator=(State const& other) {
	assert(datalayout == other.datalayout);
	assert(memptr != other.memptr);
	memcpy(memptr, other.memptr, datalayout.N*sizeof(comp));
}

inline void State::zero() {
	memset(memptr, 0x00, datalayout.N*sizeof(comp));
}

inline void State::set_by_func(comp (*initfunc)(double, double)) {
	for (size_t y=0; y<datalayout.sizey; y++) {
		double dy = datalayout.get_posy(y);
		for (size_t x=0; x<datalayout.sizex; x++) {
			double dx = datalayout.get_posy(x);
			(*this)(x,y) = initfunc(dx,dy);
		}
	}
}

inline void State::normalize(double target_norm) {
	(*this) *= (target_norm / (*this).norm());
}

// Arithmetic with other States

inline State& State::operator+=(const State& other) {
	assert(datalayout == other.datalayout);
	const comp mult = 1;
	cblas_zaxpy(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(&mult),
			reinterpret_cast<const double*>(other.memptr), 1,
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

inline State& State::operator-=(const State& other) {
	assert(datalayout == other.datalayout);
	const comp mult = -1;
	cblas_zaxpy(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(&mult),
			reinterpret_cast<const double*>(other.memptr), 1,
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

// Arithmetic with arrays

inline State& State::operator*=(const double& other) {
	cblas_zdscal(static_cast<int>(datalayout.N), other,
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

inline State& State::operator/=(const double& other) {
	const double c = 1.0/other;
	cblas_zdscal(static_cast<int>(datalayout.N), c,
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

inline State& State::operator*=(const comp& other) {
	cblas_zscal(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(&other),
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

inline State& State::operator/=(const comp& other) {
	const comp c = comp(1)/other;
	cblas_zscal(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(&c),
			reinterpret_cast<double*>(this->memptr), 1);
	return *this;
}

template <typename Type>
inline void State::pointwise_multiply(Type const* values) {
	for (size_t i=0; i<datalayout.N; i++)
		memptr[i] *= values[i];
}

template <typename Type>
inline void State::pointwise_divide(Type const* values) {
	for (size_t i=0; i<datalayout.N; i++)
		memptr[i] /= values[i];
}

// Multiply by a purely imaginary array (the imaginary part is given by argument 'values'), shiting
// the x-coordinates by one in the multiplication.
template <typename Type>
inline void State::pointwise_multiply_imaginary_shiftx(Type const* values) {
	for (size_t x=datalayout.sizex-1; x>0; x--)
		for (size_t y=0; y<datalayout.sizey; y++)
			datalayout.value(memptr, x, y) = datalayout.value(memptr, x-1, y)*comp(0, datalayout.value(values, x-1, y));
	for (size_t y=0; y<datalayout.sizey; y++)
		datalayout.value(memptr, 0, y) = 0;
}

// Multiply each constant-y-coordinate slice with a y-dependent value
template <typename Type>
inline void State::pointwise_multiply_y(Type const* values) {
	for (size_t y=0; y<datalayout.sizey; y++) {
		const Type val = values[y];
		for (size_t x=0; x<datalayout.sizex; x++)
			datalayout.value(memptr, x, y) *= val;
	}
}

template <typename Type>
inline void State::pointwise_multiply_and_add(Type const* values, const State& addstate) {
	assert(datalayout == addstate.datalayout);
	for (size_t i=0; i<datalayout.N; i++)
		memptr[i] = memptr[i]*values[i] + addstate.memptr[i];
}

// Dot product and norm

inline double State::norm() const {
	return cblas_dznrm2(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(this->memptr), 1)*datalayout.dx;
}

inline comp State::dot(const State& other) const {
	assert(datalayout == other.datalayout);
	comp sum;
	cblas_zdotc_sub(static_cast<int>(datalayout.N),
			reinterpret_cast<const double*>(this->memptr), 1,
			reinterpret_cast<const double*>(other.memptr), 1,
			reinterpret_cast<double _Complex*>(&sum));
	return sum*datalayout.dx*datalayout.dx;
}

// The free functions that implement basic arithmetic

// With states
inline const State operator+(const State& lhand, const State& rhand) {
	return State(lhand) += rhand;
}

inline const State operator-(const State& lhand, const State& rhand) {
	return State(lhand) -= rhand;
}

// With other types

template<typename Type>
inline const State operator*(const State& lhand, const Type& rhand) {
	return State(lhand) *= rhand;
}

template<typename Type>
inline const State operator/(const State& lhand, const Type& rhand) {
	return State(lhand) /= rhand;
}

template<typename Type>
inline const State operator*(const Type& lhand, const State& rhand) {
	return State(rhand) *= lhand;
}

template<typename Type>
inline const State operator/(const Type& lhand, const State& rhand) {
	return State(rhand) /= lhand;
}

#endif // _STATES_HPP_
