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

#include "test_parser.hpp"

using namespace std;

TEST(parser, parser) {
	vector<double> empty;
	double testvalues_array[4] = {1.0, 2.0, 3.0, 4.0};
	vector<double> testvalues1(testvalues_array, testvalues_array+1);
	vector<double> testvalues2(testvalues_array, testvalues_array+2);
	vector<double> testvalues4(testvalues_array, testvalues_array+4);
	// Test that empty string throws
	EXPECT_THROW(parse_parameter_string(""), ParseError);
	// Test correct whitespace stripping
	EXPECT_EQ(make_pair(string("foo"),empty) , parse_parameter_string("foo"));
	EXPECT_EQ(make_pair(string("bar"),empty) , parse_parameter_string("    bar"));
	EXPECT_EQ(make_pair(string("baz"),empty) , parse_parameter_string("baz   "));
	EXPECT_EQ(make_pair(string("foo"),empty) , parse_parameter_string("foo()"));
	EXPECT_EQ(make_pair(string("foo"),empty) , parse_parameter_string(" \t   foo()   "));
	EXPECT_EQ(make_pair(string("foo"),empty) , parse_parameter_string("    foo() \n  "));
	// Tests with nonempty values
	EXPECT_EQ(make_pair(string("foo"),testvalues1) , parse_parameter_string("foo(1.0)"));
	EXPECT_EQ(make_pair(string("frobs"),testvalues2) , parse_parameter_string("frobs(1.0, 2.0)"));
	EXPECT_EQ(make_pair(string("bar"),testvalues4) , parse_parameter_string(" bar(1.0,2.0,   3.0,\t4.0)\n"));
	// Some more throw tests
	EXPECT_THROW(parse_parameter_string("foo("), ParseError);
	EXPECT_THROW(parse_parameter_string("foo)("), ParseError);
	EXPECT_THROW(parse_parameter_string("foo bar"), ParseError);
	EXPECT_THROW(parse_parameter_string("foo(bar"), ParseError);
	EXPECT_THROW(parse_parameter_string("    foo  (1.0)   "), ParseError);
	EXPECT_THROW(parse_parameter_string("foo(1.0, 2.0) --baz"), ParseError);
	EXPECT_THROW(parse_parameter_string("foo(1.0, 2.0,)"), ParseError);
}
