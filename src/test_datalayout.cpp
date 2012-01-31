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
 * Unit tests for the datalayout class.
 */

#include "test_datalayout.hpp"

TEST(datalayout, construction) {
	const DataLayout dl(16, 8, 0.5);
	ASSERT_EQ(dl.sizex, 16u);
	ASSERT_EQ(dl.sizey, 8u);
	ASSERT_EQ(dl.dx, 0.5);
	ASSERT_EQ(dl.lenx, 8.0);
	ASSERT_EQ(dl.leny, 4.0);
}

TEST(datalayout, position_values) {
	const DataLayout dl(2, 2, 0.5);
	const DataLayout dl2(5, 3, 0.5);
	ASSERT_EQ(dl.get_posx(0), -0.25);
	ASSERT_EQ(dl.get_posx(1), 0.25);
	ASSERT_EQ(dl.get_posy(0), -0.25);
	ASSERT_EQ(dl.get_posy(1), 0.25);
	ASSERT_EQ(dl2.get_posx(0), -1.0);
	ASSERT_EQ(dl2.get_posx(1), -0.5);
	ASSERT_EQ(dl2.get_posx(2), 0);
	ASSERT_EQ(dl2.get_posx(3), 0.5);
	ASSERT_EQ(dl2.get_posx(4), 1.0);
	ASSERT_EQ(dl2.get_posy(0), -0.5);
	ASSERT_EQ(dl2.get_posy(1), 0);
	ASSERT_EQ(dl2.get_posy(2), 0.5);
}

TEST(datalayout, values_method) {
	const DataLayout dl(2, 2, 1.0);
	double* a = new double[4];
	for (size_t n=0; n<4; n++) {
		a[n] = static_cast<double>(n);
	}
	ASSERT_TRUE(typeid(dl.value(a, 0, 0)) == typeid(a[0])); 
	ASSERT_EQ(dl.value(a, 0, 0), 0);
	ASSERT_EQ(dl.value(a, 1, 0), 1);
	ASSERT_EQ(dl.value(a, 0, 1), 2);
	ASSERT_EQ(dl.value(a, 1, 1), 3);
	double const* ac = a;
	ASSERT_TRUE(typeid(dl.value(ac, 0, 0)) == typeid(a[0])); 
	ASSERT_EQ(dl.value(ac, 0, 0), 0);
	ASSERT_EQ(dl.value(ac, 1, 0), 1);
	ASSERT_EQ(dl.value(ac, 0, 1), 2);
	ASSERT_EQ(dl.value(ac, 1, 1), 3);
	delete[] a;
}
