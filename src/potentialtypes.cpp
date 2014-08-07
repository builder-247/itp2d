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

const double HarmonicOscillator::default_frequency = 1;
const double HarmonicOscillator::default_x0 = 0;
const double HarmonicOscillator::default_y0 = 0;
const double EllipticOscillator::default_frequency_x = 1;
const double EllipticOscillator::default_frequency_y = 1.6180339887498948482;
const double PrettyHardSquare::default_exponent = 8;
const double HenonHeiles::default_a = 205.0/42.0;
const double HenonHeiles::default_b = -13.0/3.0;
const double GaussianPotential::default_amplitude = 1;
const double GaussianPotential::default_width = 1;
const double GaussianPotential::default_x0 = 0;
const double GaussianPotential::default_y0 = 0;
const double QuarticPotential::default_b = 0.01;
const double SquareOscillator::default_alpha = 8;
const double PowerOscillator::default_exponent = 4;
const double PowerOscillator::default_w = 1;
const double RingPotential::default_radius = 3;
const double RingPotential::default_width = 1;
const double RingPotential::default_exponent = 2;
const double RingPotential::default_asymm_amplitude = 0;
const double RingPotential::default_asymm_width = 1;
const double CoshPotential::default_amplitude = 1;
const double CoshPotential::default_length_scale = 1;
const double SoftStadium::default_radius = 1;
const double SoftStadium::default_center_length = 2;
const double SoftStadium::default_height = 100;
const double SoftStadium::default_a = 1;
const double SoftStadium::default_b = 10;
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
		return new HarmonicOscillator(params);
	if (name == "elliptic")
		return new EllipticOscillator(params);
	if (name == "prettyhardsquare")
		return new PrettyHardSquare(params);
	if (name == "softpentagon")
		return new SoftPentagon(params);
	if (name == "henonheiles" or name == "henon")
		return new HenonHeiles(params);
	if (name == "gaussian" or name == "gaussianblob")
		return new GaussianPotential(params);
	if (name == "quartic" or name == "quarticoscillator")
		return new QuarticPotential(params);
	if (name == "squareoscillator")
		return new SquareOscillator(params);
	if (name == "power" or name == "poweroscillator")
		return new PowerOscillator(params);
	if (name == "ring" or name == "ringpotential")
		return new RingPotential(params);
	if (name == "cosh" or name == "coshpotential")
		return new CoshPotential(params);
	if (name == "softstadium")
		return new SoftStadium(params);
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
HarmonicOscillator::HarmonicOscillator(double omega, double orig_x, double orig_y) :
		w(omega),
		x0(orig_x),
		y0(orig_y) {
	init();
}

HarmonicOscillator::HarmonicOscillator(std::vector<double> params) {
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

void HarmonicOscillator::init() {
	if (w < 0)
		throw InvalidPotentialType("harmonic oscillator with negative frequency");
	std::stringstream ss;
	ss << "harmonic(" << w << ")";
	description = ss.str();
}

// EllipticOscillator

EllipticOscillator::EllipticOscillator(double omega_x, double omega_y) :
		wx(omega_x),
		wy(omega_y) {
	init();
}

EllipticOscillator::EllipticOscillator(std::vector<double> params) {
	if (params.empty()) {
		wx = default_frequency_x;
		wy = default_frequency_y;
	}
	else if (params.size() == 2) {
		wx = params[0];
		wy = params[1];
	}
	else
		throw InvalidPotentialType("elliptic oscillator potential takes either zero or three parameters");
	init();
}

void EllipticOscillator::init() {
	if (wx < 0 or wy < 0)
		throw InvalidPotentialType("elliptic oscillator with negative frequency");
	std::stringstream ss;
	ss << "elliptic(" << wx <<", " << wy << ")";
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

// QuarticPotential

QuarticPotential::QuarticPotential(std::vector<double> params) {
	if (params.empty()) {
		b = default_b;
	}
	else if (params.size() == 1) {
		b = params[0];
	}
	else
		throw InvalidPotentialType("quartic potential only takes one parameter");
	init();
}

void QuarticPotential::init() {
	if (b < 0)
		throw InvalidPotentialType("quartic potential with negative parameter");
	std::stringstream ss;
	ss << "quartic(" << b << ")";
	description = ss.str();
}

// SquareOscillator

SquareOscillator::SquareOscillator(std::vector<double> params) {
	if (params.empty())
		alpha = default_alpha;
	else if (params.size() == 1)
		alpha = params[0];
	else
		throw InvalidPotentialType("square oscillator potential only takes one parameter");
	init();
}

void SquareOscillator::init() {
	std::stringstream ss;
	ss << "squareoscillator(" << alpha << ")";
	description = ss.str();
}

// PowerOscillator

PowerOscillator::PowerOscillator(double exponent, double omega) :
		a(exponent),
		w(omega) {
	init();
}

PowerOscillator::PowerOscillator(std::vector<double> params) {
	if (params.empty()) {
		a = default_exponent;
		w = default_w;
	}
	else if (params.size() == 1) {
		a = params[0];
		w = default_w;
	}
	else if (params.size() == 2) {
		a = params[0];
		w = params[1];
	}
	else
		throw InvalidPotentialType("power oscillator potential takes either zero, one or two parameters");
	init();
}

void PowerOscillator::init() {
	if (w < 0)
		throw InvalidPotentialType("power oscillator with negative \"frequency\"");
	std::stringstream ss;
	ss << "poweroscillator(" << a << "," << w << ")";
	description = ss.str();
}

// RingPotential

RingPotential::RingPotential(double radius, double width, double exponent, double asymm_amplitude, double asymm_width) :
		r(radius),
		w(width),
		e(exponent),
		asymm_A(asymm_amplitude),
		asymm_w(asymm_width) {
	init();
}

RingPotential::RingPotential(std::vector<double> params) {
	if (params.empty()) {
		r = default_radius;
		w = default_width;
		e = default_exponent;
		asymm_A = default_asymm_amplitude;
		asymm_w = default_asymm_width;
	}
	else if (params.size() == 3) {
		r = params[0];
		w = params[1];
		e = params[2];
		asymm_A = default_asymm_amplitude;
		asymm_w = default_asymm_width;
	}
	else if (params.size() == 5) {
		r = params[0];
		w = params[1];
		e = params[2];
		asymm_A = params[3];
		asymm_w = params[4];
	}
	else
		throw InvalidPotentialType("ring potential takes either zero, three, or five parameters");
	init();
}

void RingPotential::init() {
	if (r < 0)
		throw InvalidPotentialType("ring oscillator with negative radius");
	if (w <= 0)
		throw InvalidPotentialType("ring oscillator with non-positive width");
	std::stringstream ss;
	ss << "ring(" << r << "," << w << "," << e << "," << asymm_A << "," << asymm_w << ")";
	description = ss.str();
}

// CoshPotential

CoshPotential::CoshPotential(double amplitude, double length_scale) : A(amplitude), L(length_scale) {
	init();
}

CoshPotential::CoshPotential(std::vector<double> params) {
	if (params.empty()) {
		A = default_amplitude;
		L = default_length_scale;
	}
	else if (params.size() == 2) {
		A = params[0];
		L = params[1];
	}
	else
		throw InvalidPotentialType("cosh potential takes either zero or two parameters");
	init();
}

void CoshPotential::init() {
	if (A <= 0)
		throw InvalidPotentialType("cosh potential with non-positive amplitude");
	if (L <= 0)
		throw InvalidPotentialType("cosh potential with non-positive length scale");
	std::stringstream ss;
	ss << "cosh(" << A << "," << L << ")";
	description = ss.str();
}

// SoftStadium

SoftStadium::SoftStadium(double radius, double center_length, double height, double _a, double _b) :
		R(radius), halfL(center_length/2), V(height), a(_a), b(_b) {
	init();
}

SoftStadium::SoftStadium(std::vector<double> params) {
	if (params.empty()) {
		R = default_radius;
		halfL = default_center_length/2;
		V = default_height;
		a = default_a;
		b = default_b;
	}
	else if (params.size() == 3) {
		R = params[0];
		halfL = params[1];
		V = params[2];
		a = default_a;
		b = default_b;
	}
	else if (params.size() == 5) {
		R = params[0];
		halfL = params[1];
		V = params[2];
		a = params[3];
		b = params[4];
	}
	else
		throw InvalidPotentialType("soft stadium potential takes either zero, three or five parameters");
	init();
}

void SoftStadium::init() {
	if (R <= 0)
		throw InvalidPotentialType("soft stadium with non-positive radius");
	if (halfL <= 0)
		throw InvalidPotentialType("soft stadium with non-positive center length");
	if (V <= 0)
		throw InvalidPotentialType("soft stadium with non-positive height");
	if (a <= 0)
		throw InvalidPotentialType("soft stadium with non-positive a parameter");
	if (b <= 0)
		throw InvalidPotentialType("soft stadium with non-positive b parameter");
	std::stringstream ss;
	ss << "softstadium(" << R << "," << halfL*2 << "," << V << "," << a << "," << b << ")";
	description = ss.str();
}
