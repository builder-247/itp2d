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
 * Unit tests for the simple operator classes. This mainly tests whether the
 * implementations of OperatorProduct and OperatorSum work correctly.
 */

#include "test_operators.hpp"

class operators : public testing::Test {
public:
	operators() : dl(2, 2, 1.0), work(2, dl), A(dl), OriginalA(dl) {
		A(0,0) = comp(0,0);
		A(0,1) = comp(0,1);
		A(1,0) = comp(1,0);
		A(1,1) = comp(1,1);
		OriginalA = A;
	}
	const DataLayout dl;
	StateArray work;
	State A;
	State OriginalA;
};

TEST_F(operators, nulloperator) {
	NullOperator I;
	I(A, work);
	EXPECT_EQ(A, OriginalA);
}

TEST_F(operators, zerooperator) {
	ZeroOperator O;
	O(A, work);
	State Zero(dl);
	Zero.zero();
	EXPECT_EQ(A, Zero);
}

TEST_F(operators, constantoperator) {
	ConstantOperator Two(2);
	Two(A, work);
	EXPECT_EQ(A, 2*OriginalA);
	A = OriginalA;
	ConstantOperator I(comp(0,1));
	I(A, work);
	EXPECT_EQ(A, comp(0,1)*OriginalA);
}

TEST_F(operators, operatorproduct) {
	ConstantOperator Two(2);
	ConstantOperator I(comp(0,1));
	OperatorProduct Prod = Two*Two*I;
	Prod(A, work);
	EXPECT_EQ(A, comp(0,4)*OriginalA);
}

TEST_F(operators, operatorsum) {
	ConstantOperator Two(2);
	ConstantOperator I(comp(0,1));
	OperatorSum Sum = Two+Two+I;
	Sum(A, work);
	EXPECT_EQ(A, 2*OriginalA + 2*OriginalA + comp(0,1)*OriginalA);
}

TEST_F(operators, matrixelement) {
	ConstantOperator Two(2);
	ConstantOperator I(comp(0,1));
	EXPECT_EQ(Two.matrixelement(A, A, work)/A.dot(A), comp(2));
	EXPECT_EQ(I.matrixelement(A, A, work)/A.dot(A), comp(0,1));
}
