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
 * Unit tests for the command line parser. This is basically a stub. More tests
 * should be implemented.
 */

#include "test_commandlineparser.hpp"

class commandlineparser : public testing::Test {
	public:
	CommandLineParser parser;
};

TEST_F(commandlineparser, highmem) {
	std::vector<std::string> fakeargv(2);
	fakeargv[0] = "test";
	fakeargv[1] = "--highmem";
	parser.parse(fakeargv);
	ASSERT_EQ(parser.get_params().get_ortho_algorithm(), HighMem);
}

// TODO: Add unit tests to other features of the command line parser
