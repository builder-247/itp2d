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

#include "parser.hpp"
using namespace std;

// If you make changes, remember to also update the in-line documentation in
// commandlineparser.cpp


pair<string,vector<double> > parse_parameter_string(string const& s) {
	// First, search for whitespace around the string so we can ignore it
	const size_t start = s.find_first_not_of(whitespace_characters);
	if (start == string::npos)
		throw ParseError(s); // string was completely whitespace
	const size_t end = s.find_last_not_of(whitespace_characters);
	// Then, find the first opening parentheses or whitespace. Everything
	// between the first non-whitespace and this point is considered to be the
	// name
	const size_t name_end = s.find_first_of(name_delim_characters, start);
	string name(s, start, name_end-start);
	// OK, then the values
	vector<double> values;
	// Check that the rest of the string is enclosed in parentheses, or completely
	// whitespace
	if (s.find_first_not_of(whitespace_characters, name_end) == string::npos)
		return make_pair(name, values); // was completely whitespace
	if (not (s.at(name_end) == '(' and s.at(end) == ')'))
		throw ParseError(s);	// was not enclosed in parentheses
	// Parse values; first create a istringstream with the part enclosed in the
	// parentheses
	string valuestring(s, name_end+1, end-name_end-1);
	if (valuestring.empty())
		return make_pair(name,values);
	istringstream in(valuestring);
	double val;
	while (not in.eof()) {
		in >> val;
		if (in.fail())
			throw ParseError(s);
		values.push_back(val);
		in.ignore(32, ',');
	}
	return make_pair(name, values);
}
