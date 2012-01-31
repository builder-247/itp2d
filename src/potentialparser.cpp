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

#include "potentialparser.hpp"

// If you make changes, remember to also update the in-line documentation in
// commandlineparser.cpp

PotentialType const* parse_potential_description(std::string const& str) {
	/* This could be (and has been) done with regular expressions, but since
	 * this functionality is not commonly available on standard C++, it seems
	 * silly to include a huge dependency for a very minor part of the
	 * program.
	 */
	// simple cases first
	if (str == "zero")
		return new ZeroPotential();
	if (str == "harmonic" or str == "default")
		return new HarmonicPotential();
	if (str == "prettyhardsquare")
		return new PrettyHardSquare();
	if (str == "softpentagon")
		return new SoftPentagon();
	if (str == "henonheiles")
		return new HenonHeiles();
	if (str == "gaussian")
		return new Gaussian();
	// then potentials with parameters
	/* For some reason this is simpler to code in C than C++. */
	double param, param2, param3, param4;
	if (sscanf(str.c_str(), "harmonic(%lf)", &param))
		return new HarmonicPotential(param);
	if (sscanf(str.c_str(), "prettyhardsquare(%lf)", &param))
		return new PrettyHardSquare(param);
	if (sscanf(str.c_str(), "henonheiles(%lf,%lf)", &param, &param2))
		return new HenonHeiles(param, param2);
	if (sscanf(str.c_str(), "gaussian(%lf,%lf)", &param, &param2))
		return new Gaussian(param, param2);
	if (sscanf(str.c_str(), "gaussian(%lf,%lf,%lf,%lf)", &param, &param2, &param3, &param4))
		return new Gaussian(param, param2, param3, param3);
	// if nothing matches, throw
	throw UnknownPotentialType(str);
}
