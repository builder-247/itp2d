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

/* Potential type parser that returns a PotentialType* from an user supplied
 * string. Remember to delete the pointer you get.
 */

#ifndef _POTENTIALPARSER_HPP_
#define _POTENTIALPARSER_HPP_

#include <string>

#include "potentialtypes.hpp"
#include "exceptions.hpp"

PotentialType const* parse_potential_description(std::string const& str);

#endif // _POTENTIALPARSER_HPP_
