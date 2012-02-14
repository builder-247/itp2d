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
 * Unit tests for the StateArray class.
 */

#include "test_statearray.hpp"

class statearray : public testing::Test {
public:
	statearray() : dl(2,2,1.0) {
		dataptr = reinterpret_cast<comp*>(fftw_malloc(4*3*sizeof(comp)));
		for (int i=0; i<4*3; i++ )
			dataptr[i] = comp(i,pi);
	}
	~statearray() {
		fftw_free(dataptr);
	}
	const DataLayout dl;
	comp* dataptr;
};

TEST_F(statearray, construction) {
	StateArray A(3, dl);
	EXPECT_TRUE(&A.datalayout == &dl);
	EXPECT_TRUE(A.get_dataptr() != NULL);
	EXPECT_EQ(A.size(), 3u);
	StateArray B(3, dl, dataptr);
	EXPECT_EQ(B.size(), 3u);
	EXPECT_TRUE(&B.datalayout == &dl);
	EXPECT_EQ(B.get_dataptr(), dataptr);
	EXPECT_EQ(B[0](0,0), comp(0,pi));
	EXPECT_EQ(B[0](1,0), comp(1,pi));
	EXPECT_EQ(B[0](0,1), comp(2,pi));
	EXPECT_EQ(B[0](1,1), comp(3,pi));
	EXPECT_EQ(B[1](0,0), comp(4,pi));
	EXPECT_EQ(B[1](1,0), comp(5,pi));
	EXPECT_EQ(B[1](0,1), comp(6,pi));
	EXPECT_EQ(B[1](1,1), comp(7,pi));
	EXPECT_EQ(B[2](0,0), comp(8,pi));
	EXPECT_EQ(B[2](1,0), comp(9,pi));
	EXPECT_EQ(B[2](0,1), comp(10,pi));
	EXPECT_EQ(B[2](1,1), comp(11,pi));
}

TEST_F(statearray, slicing) {
	StateArray A(3, dl, dataptr);
	StateArray S0(A,0);
	StateArray S1(A,1);
	StateArray S2(A,2);
	StateArray S3(A,3);
	EXPECT_EQ(S0.size(), A.size());
	EXPECT_EQ(S1.size(), A.size()-1);
	EXPECT_EQ(S2.size(), A.size()-2);
	EXPECT_EQ(S3.size(), A.size()-3);
	EXPECT_TRUE(A[0] != A[1]);
	EXPECT_EQ(A[0].data_ptr(), S0[0].data_ptr());
	EXPECT_EQ(A[1].data_ptr(), S0[1].data_ptr());
	EXPECT_EQ(A[2].data_ptr(), S0[2].data_ptr());
	EXPECT_EQ(A[1].data_ptr(), S1[0].data_ptr());
	EXPECT_EQ(A[2].data_ptr(), S1[1].data_ptr());
	EXPECT_EQ(A[2].data_ptr(), S2[0].data_ptr());
	StateArray S11(A,1,1);
	EXPECT_EQ(S11.size(), 1u);
	EXPECT_EQ(A[1].data_ptr(), S11[0].data_ptr());
}
