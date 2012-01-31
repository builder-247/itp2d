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

#ifndef _CONVERGENCEPARSER_HPP_
#define _CONVERGENCEPARSER_HPP_

/* Parser that returns a pointer to a ConvergenceTest instance based on a
 * description string provided by the user. Remember to delete the pointer you
 * get!
 */

#include <string>

#include "convergence.hpp"
#include "exceptions.hpp"

ConvergenceTest const* parse_convergence_description(std::string const& str);

#endif // _CONVERGENCEPARSER_HPP_
