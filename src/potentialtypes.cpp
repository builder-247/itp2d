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

#include "potentialtypes.hpp"

// Default values for parameters

const double HarmonicPotential::default_frequency = 1;
const double PrettyHardSquare::default_exponent = 8;
const double HenonHeiles::default_a = 205.0/42.0;
const double HenonHeiles::default_b = -13.0/3.0;
const double Gaussian::default_amplitude = 1;
const double Gaussian::default_width = 1;
const double Gaussian::default_x0 = 0;
const double Gaussian::default_y0 = 0;

// The parser delegator

PotentialType const* parse_potential_description(std::string const& str) {
	name_parameters_pair p;
	try {
		p = parse_parameter_string(str);
	}
	catch (ParseError& e) {
		std::cerr << e.what() << std::endl;
		throw InvalidPotentialType(str);
	}
	std::string const& name = p.first;
	std::vector<double> const& params = p.second;
	// Simply delegate to the individual constructors based on name
	// Unfortunately we cannot use switch-case here, since we are dealing with
	// strings
	if (name == "zero")
		return new ZeroPotential(params);
	if (name == "harmonic")
		return new HarmonicPotential(params);
	if (name == "prettyhardsquare")
		return new PrettyHardSquare(params);
	if (name == "softpentagon")
		return new SoftPentagon(params);
	if (name == "henonheiles" or name == "henon")
		return new HenonHeiles(params);
	if (name == "gaussian" or name == "gaussianblob")
		return new Gaussian(params);
	else
		throw UnknownPotentialType(str);
	return NULL;
}
