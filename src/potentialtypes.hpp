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

/*
 * Different potential types are defined here.
 *
 * This is where you need to make changes if you want to implement new
 * potential functions. If you add a potential, remember to also update
 * potentialtypes.cpp where the parsing function is implemented. Also make
 * modifications in commandlineparser.cpp to include your potential in the
 * in-line documentation.
 */

#ifndef _POTENTIALTYPES_HPP_
#define _POTENTIALTYPES_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "itp2d_common.hpp"
#include "exceptions.hpp"
#include "parser.hpp"

// Main interface class.
// A potential type is simply:
// 	* A function giving its value on a certain point (x,y)
// 	* A descriptive string
// In addition, you need a constructor that constructs a potential from a
// vector of doubles. This is used for creating potentials from user-provided
// parameters.

class PotentialType {
	public:
		virtual double operator()(double x, double y) const = 0;
		inline std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

// A parser function for constructing potentials from a user-provided string.
// First, the generic parse_parameter_string is used to parse a string of the
// form name(param1,param2,...) into a name and a vector of doubles. Then the
// construction is delegated to the respective subclass constructors based on
// the name.
//
// This parser could be (and has been) put into a separate file. However, it is
// good to have *everything* relating to a single potential close to each
// other.

PotentialType const* parse_potential_description(std::string const& str);

// Potential types

class ZeroPotential : public PotentialType {
	public:
		ZeroPotential() {
			init();
		}
		ZeroPotential(std::vector<double> params) {
			if (not params.empty())
				throw InvalidPotentialType("zero potential does not take parameters");
			init();
		}
		inline double operator()(__attribute__((unused)) double x, __attribute__((unused)) double y) const { return 0; }
	private:
		void init() {
			description = "zero";
		}
};

// This potential type is not used at the moment, but it is implemented here in
// case e.g. some method of reading the potential from a file is made available
// later
class UserSetPotential : public PotentialType {
	public:
		typedef double (*potfunc)(double, double);
		UserSetPotential(std::string desc, potfunc func) : f(func) {
			description = desc;
		}
		inline double operator()(double x, double y) const { return f(x,y); }
	private:
		potfunc f;
};

// The harmonic oscillator
class HarmonicPotential : public PotentialType {
	public:
		static const double default_frequency;
		static const double default_x0;
		static const double default_y0;
		HarmonicPotential(double omega=default_frequency) : w(omega) {
			init();
		}
		HarmonicPotential(std::vector<double> params) {
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
		inline double operator()(double x, double y) const {
			const double px = x-x0;
			const double py = y-y0;
			return 0.5*w*(px*px+py*py);
		}
	private:
		double w;
		double x0;
		double y0;
		void init() {
			if (w < 0)
				throw InvalidPotentialType("harmonic oscillator with negative frequency");
			std::stringstream ss;
			ss << "harmonic(" << w << ")";
			description = ss.str();
		}
};

// A pretty hard square of length pi
class PrettyHardSquare : public PotentialType {
	public:
		static const double default_exponent;
		PrettyHardSquare(double e=default_exponent) : exponent(e) {
			init();
		}
		PrettyHardSquare(std::vector<double> params) {
			if (params.empty())
				exponent = default_exponent;
			else if (params.size() == 1)
				exponent = params[0];
			else
				throw InvalidPotentialType("pretty hard square potential only takes one parameter");
			init();
		}
		inline double operator()(double x, double y) const {
			const double ax = fabs(2*x/pi);
			const double ay = fabs(2*y/pi);
			const double r = (ax > ay)? ax : ay;
			return pow(r,exponent);
		}
	private:
		double exponent;
		void init() {
			std::stringstream ss;
			ss << "prettyhardsquare(" << exponent << ")";
			description = ss.str();
		}
};

// A soft potential with a pentagon shape
class SoftPentagon : public PotentialType {
	public:
		SoftPentagon() {
			init();
		}
		SoftPentagon(std::vector<double> params) {
			if (not params.empty())
				throw InvalidPotentialType("soft pentagon potential does not take parameters");
			init();
		}
		inline double operator()(double x, double y) const {
			const double r2 = x*x+y*y;
			const double t = atan2(y,x);
			return 0.5*r2*(1 + 0.5*sin(5*t));
		}
	private:
		void init() {
			description = "softpentagon";
		}
};

// A Hénon-Heiles type potential, modified so that all trajectories are bounded and generalized
// as in http://mathworld.wolfram.com/Henon-HeilesEquation.html
class HenonHeiles : public PotentialType {
	public:
		static const double default_a;
		static const double default_b;
		HenonHeiles(double _a=default_a, double _b=default_b) : a(_a), b(_b) {
			init();
		}
		HenonHeiles(std::vector<double> params) {
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
		inline double operator()(double x, double y) const {
			const double r2 = x*x+y*y;
			const double r3 = pow(r2, 3.0/2.0);
			const double r4 = r2*r2;
			const double t = atan2(y,x);
			return r4 + a*r2 + b*r3*cos(3*t);
		}
	private:
		double a;
		double b;
		void init() {
			std::stringstream ss;
			ss << "henonheiles(" << a << "," << b << ")";
			description = ss.str();
		}
};

// A Gaussian blob with a prescribed amplitude and width, centered at (x0,y0)
class Gaussian : public PotentialType {
	public:
		static const double default_amplitude;
		static const double default_width;
		static const double default_x0;
		static const double default_y0;
		Gaussian(double _amplitude=default_amplitude, double _width=default_width,
				double _x0=default_x0, double _y0=default_y0) :
				amplitude(_amplitude), width(_width), x0(_x0), y0(_y0) {
			init();
		}
		Gaussian(std::vector<double> params) {
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
		inline double operator()(double x, double y) const {
			const double xp = x-x0;
			const double yp = y-y0;
			const double exponent = -(xp*xp + yp*yp)/(2*width*width);
			return amplitude*exp(exponent);
		}
	private:
		double amplitude;
		double width;
		double x0;
		double y0;
		void init() {
			if (amplitude < 0)
				throw InvalidPotentialType("gaussian potential with negative amplitude");
			if (width == 0)
				throw InvalidPotentialType("gaussian potential with zero width");
			std::stringstream ss;
			ss << "gaussian(" << amplitude << "," << width << "," << x0 << "," << y0 << ")";
			description = ss.str();
		}
};

#endif // _POTENTIALTYPES_HPP_
