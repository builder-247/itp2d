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

/* A generic parser function used to parse descriptions of potential, noise and
 * convergence criteria.
 *
 * Reads a string of the format "name(value1,value2,value3)", with any number
 * of values castable to double, and returns a pair of a string and a vector of
 * doubles.
 */

#ifndef _PARSER_HPP_
#define _PARSER_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include "exceptions.hpp"

static const std::string whitespace_characters = " \f\n\r\t\v";
static const std::string name_delim_characters = " \f\n\r\t\v(";

typedef std::pair<std::string,std::vector<double> > name_parameters_pair;

name_parameters_pair parse_parameter_string(std::string const& s);

#endif // _PARSER_HPP_
