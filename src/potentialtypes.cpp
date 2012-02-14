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
const double HarmonicPotential::default_x0 = 0;
const double HarmonicPotential::default_y0 = 0;
const double PrettyHardSquare::default_exponent = 8;
const double HenonHeiles::default_a = 205.0/42.0;
const double HenonHeiles::default_b = -13.0/3.0;
const double GaussianPotential::default_amplitude = 1;
const double GaussianPotential::default_width = 1;
const double GaussianPotential::default_x0 = 0;
const double GaussianPotential::default_y0 = 0;

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
		return new GaussianPotential(params);
	else
		throw UnknownPotentialType(str);
	return NULL;
}

// Zero potential

ZeroPotential::ZeroPotential(std::vector<double> params) {
	if (not params.empty())
		throw InvalidPotentialType("zero potential does not take parameters");
	init();
}

// The harmonic potential

HarmonicPotential::HarmonicPotential(std::vector<double> params) {
	if (params.empty()) {
		w = default_frequency;
		x0 = default_x0;
		y0 = default_y0;
	}
	else if (params.size() == 1) {
		w = params[0];
		x0 = default_x0;
		y0 = default_y0;
	}
	else if (params.size() == 3) {
		w = params[0];
		x0 = params[1];
		y0 = params[2];
	}
	else
		throw InvalidPotentialType("harmonic oscillator potential takes either zero, one or three parameters");
	init();
}

void HarmonicPotential::init() {
	if (w < 0)
		throw InvalidPotentialType("harmonic oscillator with negative frequency");
	std::stringstream ss;
	ss << "harmonic(" << w << ")";
	description = ss.str();
}

// PrettyHardSquare

PrettyHardSquare::PrettyHardSquare(std::vector<double> params) {
	if (params.empty())
		exponent = default_exponent;
	else if (params.size() == 1)
		exponent = params[0];
	else
		throw InvalidPotentialType("pretty hard square potential only takes one parameter");
	init();
}

void PrettyHardSquare::init() {
	std::stringstream ss;
	ss << "prettyhardsquare(" << exponent << ")";
	description = ss.str();
}

// SoftPentagon

SoftPentagon::SoftPentagon(std::vector<double> params) {
	if (not params.empty())
		throw InvalidPotentialType("soft pentagon potential does not take parameters");
	init();
}

// HenonHeiles

HenonHeiles::HenonHeiles(std::vector<double> params) {
	if (params.empty()) {
		a = default_a;
		b = default_b;
	}
	else if (params.size() == 2) {
		a = params[0];
		b = params[1];
	}
	else
		throw InvalidPotentialType("Henon Heiles potential takes either two parameters or none");
	init();
}

void HenonHeiles::init() {
	std::stringstream ss;
	ss << "henonheiles(" << a << "," << b << ")";
	description = ss.str();
}

// GaussianPotential

GaussianPotential::GaussianPotential(std::vector<double> params) {
	if (params.empty()) {
		amplitude = default_amplitude;
		width = default_width;
		x0 = default_x0;
		y0 = default_y0;
	}
	else if (params.size() == 2) {
		amplitude = params[0];
		width = params[1];
		x0 = default_x0;
		y0 = default_y0;
	}
	else if (params.size() == 4) {
		amplitude = params[0];
		width = params[1];
		x0 = params[2];
		y0 = params[3];
	}
	else
		throw InvalidPotentialType("gaussian potential takes 0, 2 or 4 parameters");
	init();
}

void GaussianPotential::init() {
	if (amplitude < 0)
		throw InvalidPotentialType("gaussian potential with negative amplitude");
	if (width == 0)
		throw InvalidPotentialType("gaussian potential with zero width");
	std::stringstream ss;
	ss << "gaussian(" << amplitude << "," << width << "," << x0 << "," << y0 << ")";
	description = ss.str();
}

