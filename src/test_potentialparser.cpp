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
 * Unit tests for the potential parser function.
 */

#include "test_potentialparser.hpp"

TEST(potentialparser, potentialparser_harmonic) {
	PotentialType const* p= parse_potential_description("harmonic");
	ASSERT_EQ(p->get_description(), "harmonic(1)");
	ASSERT_EQ((*p)(0,0), 0.0);
	ASSERT_EQ((*p)(1,1), 1.0);
	delete p;
}

TEST(potentialparser, potentialparser_harmonic_05) {
	PotentialType const* p= parse_potential_description("harmonic(0.5)");
	ASSERT_EQ(p->get_description(), "harmonic(0.5)");
	ASSERT_EQ((*p)(0,0), 0.0);
	ASSERT_EQ((*p)(1,1), 0.5);
	delete p;
}

TEST(potentialparser, potentialparser_zero) {
	PotentialType const* p= parse_potential_description("zero");
	ASSERT_EQ(p->get_description(), "zero");
	delete p;
}

TEST(potentialparser, potentialparser_prettyhardsquare) {
	PotentialType const* p= parse_potential_description("prettyhardsquare");
	ASSERT_EQ(p->get_description(), "prettyhardsquare(8)");
	delete p;
}

TEST(potentialparser, potentialparser_prettyhardsquare_628) {
	PotentialType const* p= parse_potential_description("prettyhardsquare(6.28)");
	ASSERT_EQ(p->get_description(), "prettyhardsquare(6.28)");
	delete p;
}
