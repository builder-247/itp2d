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

#include "noiseparser.hpp"

// If you make changes, remember to also update the in-line documentation in
// commandlineparser.cpp

Noise const* parse_noise_description(std::string const& str) {
	/* This could be (and has been) done with regular expressions, but since
	 * this functionality is not commonly available on standard C++, it seems
	 * silly to include a huge dependency for a very insignificant part of the
	 * program.
	 */
	// Simple cases first
	if (str == "none" or str == "no" or str == "zero")
		return new NoNoise();
	// ...then noise makers with parameters
	/* For some reason this is simpler to code in C than C++. */
	double param[5];
	if (sscanf(str.c_str(), "gaussians(%lf,%lf,%lf,%lf,%lf)", &param[0], &param[1], &param[2], &param[3], &param[4]))
		return new GaussianSpikes(param[0], param[1], param[2], param[3], param[4]);
	// if nothing matches, throw
	throw UnknownNoiseType(str);
}
