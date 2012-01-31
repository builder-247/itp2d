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

#include "convergenceparser.hpp"

/* Same story as with the potential type parser: this could be implemented with
 * regular expressions, but that would incur a heavy dependancy for a simple
 * thing. */

// If you make changes, remember to also update the in-line documentation in
// commandlineparser.cpp

ConvergenceTest const* parse_convergence_description(std::string const& str) {
	// simple cases first
	if (str == "no" or str == "none" or str == "null")
		return new NoConvergenceTest();
	if (str == "onestep" or str == "one-step")
		return new OneStepConvergenceTest();
	// then potentials with parameters
	double param, param2;
	if (sscanf(str.c_str(), "absEchange(%lf)", &param) or (sscanf(str.c_str(), "absEdelta(%lf)", &param)))
		return new AbsoluteEnergyChangeTest(param);
	if (sscanf(str.c_str(), "relEchange(%lf)", &param) or (sscanf(str.c_str(), "relEdelta(%lf)", &param)))
		return new RelativeEnergyChangeTest(param);
	if (sscanf(str.c_str(), "deviation(%lf,%lf)", &param, &param2))
		return new EnergyDeviationChangeTest(param, param2);
	// if nothing matches, throw
	throw UnknownConvergenceType(str);
}
