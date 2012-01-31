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

/* Noise type parser that returns a Noise* from an user supplied
 * string. Remember to delete the pointer you get!
 */

#ifndef _NOISEPARSER_HPP_
#define _NOISEPARSER_HPP_

#include <string>

#include "noise.hpp"
#include "exceptions.hpp"

Noise const* parse_noise_description(std::string const& str);

#endif // _NOISEPARSER_HPP_
