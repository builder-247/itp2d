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
 * Low-level unit tests for the State class.
 */

#include "test_state.hpp"

namespace test_states_reference {
	comp initfunc(double x, double y) {
		return comp(x,y);
	}
}

struct states : public testing::Test {
public:
	states() : dl(2,2,1.0), A(dl), B(dl) {
		A(0,0) = comp(1,2);
		A(0,1) = comp(1,2);
		A(1,0) = comp(1,2);
		A(1,1) = comp(1,2);
		B(0,0) = comp(1,0);
		B(0,1) = comp(1,2);
		B(1,0) = comp(1,0);
		B(1,1) = comp(1,2);
	}
	const DataLayout dl;
	State A, B;
};

TEST_F(states, construction) {
	B = A;
	ASSERT_EQ(A, B);
	State C(dl, test_states_reference::initfunc);
	ASSERT_EQ(C(0,0), comp(dl.get_posx(0), dl.get_posy(0)));
	ASSERT_EQ(C(0,1), comp(dl.get_posx(0), dl.get_posy(1)));
	ASSERT_EQ(C(1,0), comp(dl.get_posx(1), dl.get_posy(0)));
	ASSERT_EQ(C(1,1), comp(dl.get_posx(1), dl.get_posy(1)));
	comp* dataptr = reinterpret_cast<comp*>(fftw_malloc(4*sizeof(comp)));
	dataptr[0] = comp(0,pi);
	dataptr[1] = comp(1,pi);
	dataptr[2] = 0;
	dataptr[3] = 0;
	State D(dl, dataptr);
	ASSERT_EQ(D(0,0), comp(0,pi));
	ASSERT_EQ(D(1,0), comp(1,pi));
	State E(dl, test_states_reference::initfunc, dataptr);
	ASSERT_EQ(dataptr[0], comp(dl.get_posx(0), dl.get_posy(0)));
	ASSERT_EQ(dataptr[1], comp(dl.get_posx(1), dl.get_posy(0)));
	ASSERT_EQ(dataptr[2], comp(dl.get_posx(0), dl.get_posy(1)));
	ASSERT_EQ(dataptr[3], comp(dl.get_posx(1), dl.get_posy(1)));
	fftw_free(dataptr);
}

TEST_F(states, comparison) {
	ASSERT_EQ(A, A);
	ASSERT_EQ(B, B);
	ASSERT_EQ((A!=B), (B!=A));
	ASSERT_NE((A==B), (A!=B));
	ASSERT_EQ((A==B), (B==A));
	State C(dl);
	C = B;
	C(0,0) = comp(0,pi);
	ASSERT_NE(C, B);
	ASSERT_EQ(rms_distance(A,A), 0);
	ASSERT_EQ(rms_distance(B,B), 0);
	ASSERT_EQ(rms_distance(A,B), sqrt(2));
	ASSERT_EQ(max_distance(A,A), 0);
	ASSERT_EQ(max_distance(B,B), 0);
	ASSERT_EQ(max_distance(A,B), 2);
}

// Test simple assign
TEST_F(states, assign) {
	ASSERT_EQ(A(0,0), comp(1,2));
	ASSERT_EQ(A(0,1), comp(1,2));
	ASSERT_EQ(A(1,0), comp(1,2));
	ASSERT_EQ(A(1,1), comp(1,2));
	ASSERT_NE(A(0,0), comp(1,1));
	State C = A;
	ASSERT_EQ(C(0,0), comp(1,2));
	ASSERT_EQ(C(0,1), comp(1,2));
	ASSERT_EQ(C(1,0), comp(1,2));
	ASSERT_EQ(C(1,1), comp(1,2));
	A(0,0) = 0;
	ASSERT_EQ(A(0,0), comp(0));
	A(0,1) = pi;
	ASSERT_EQ(A(0,1), pi);
}


// Test multiplying with a constant
TEST_F(states, multiply_with_constant) {
	B = 2*A;
	ASSERT_EQ(B(0,0), comp(2,4));
	ASSERT_EQ(B(0,1), comp(2,4));
	ASSERT_EQ(B(1,0), comp(2,4));
	ASSERT_EQ(B(1,1), comp(2,4));
	B = comp(0,1)*A;
	ASSERT_EQ(B(0,0), comp(-2,1));
	ASSERT_EQ(B(0,1), comp(-2,1));
	ASSERT_EQ(B(1,0), comp(-2,1));
	ASSERT_EQ(B(1,1), comp(-2,1));
}


// Test addition and subtraction with constants

TEST_F(states, addition) {
	State C(dl);
	C = A+B;
	ASSERT_EQ(C(0,0), comp(2,2));
	ASSERT_EQ(C(0,1), comp(2,4));
	ASSERT_EQ(C(1,0), comp(2,2));
	ASSERT_EQ(C(1,1), comp(2,4));
	C = A-B;
	ASSERT_EQ(C(0,0), comp(0,2));
	ASSERT_EQ(C(0,1), comp(0,0));
	ASSERT_EQ(C(1,0), comp(0,2));
	ASSERT_EQ(C(1,1), comp(0,0));
}

TEST_F(states, dot_product) {
	ASSERT_EQ(A.dot(B), comp(12,-4));
	ASSERT_EQ(A.dot(B), conj(B.dot(A)));
}

TEST_F(states, norm) {
	ASSERT_EQ(A.norm(), sqrt(20));
	ASSERT_EQ(B.norm(), sqrt(12));
}

TEST_F(states, pointwise_multiply) {
	double darray[4] = { 1.0, 2.0, 3.0, 4.0 };
	A.pointwise_multiply(darray);
	ASSERT_EQ(A(0,0), comp(1,2));
	ASSERT_EQ(A(0,1), comp(3,6));
	ASSERT_EQ(A(1,0), comp(2,4));
	ASSERT_EQ(A(1,1), comp(4,8));
}

TEST_F(states, pointwise_multiply_imaginary_shiftx) {
	double darray[4] = { 1.0, NaN, 2.0, NaN };
	A.pointwise_multiply_imaginary_shiftx(darray);
	ASSERT_EQ(A(0,0), comp(0));
	ASSERT_EQ(A(0,1), comp(0));
	ASSERT_EQ(A(1,0), comp(-2,1));
	ASSERT_EQ(A(1,1), comp(-4,2));
}

TEST_F(states, pointwise_multiply_y) {
	double darray[3] = { 1.0, 2.0, NaN };
	A.pointwise_multiply_y(darray);
	ASSERT_EQ(A(0,0), comp(1,2));
	ASSERT_EQ(A(0,1), comp(2,4));
	ASSERT_EQ(A(1,0), comp(1,2));
	ASSERT_EQ(A(1,1), comp(2,4));
}

TEST_F(states, pointwise_multiply_and_add) {
	double darray[4] = { 1.0, 2.0, 3.0, 4.0 };
	A.pointwise_multiply_and_add(darray, B);
	ASSERT_EQ(A(0,0), comp(2,2));
	ASSERT_EQ(A(0,1), comp(4,8));
	ASSERT_EQ(A(1,0), comp(3,4));
	ASSERT_EQ(A(1,1), comp(5,10));
}
