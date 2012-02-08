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
 * potential functions. Also make modifications in commandlineparser.cpp to
 * include your potential in the in-line documentation.
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
		HarmonicPotential(double omega=1) : w(omega) {}
		HarmonicPotential(std::vector<double> params) {
			if (params.size() == 0)
				w = 1.0;
			else if (params.size() == 1)
				w = params[0];
			else
				throw InvalidPotentialType("harmonic oscillator potential only takes one parameter");
		}
		inline double operator()(double x, double y) const { return 0.5*w*(x*x+y*y); }
	private:
		void init() {
			if (w < 0)
				throw InvalidPotentialType("harmonic oscillator with negative frequency");
			std::stringstream ss;
			ss << "harmonic(" << w << ")";
			description = ss.str();
		}
		double w;
};

// A pretty hard square of length pi
class PrettyHardSquare : public PotentialType {
	public:
		PrettyHardSquare(double e=8) : exponent(e) {
			std::stringstream ss;
			ss << "prettyhardsquare(" << e << ")";
			description = ss.str();
		}
		inline double operator()(double x, double y) const {
			const double ax = fabs(2*x/pi);
			const double ay = fabs(2*y/pi);
			const double r = (ax > ay)? ax : ay;
			return pow(r,exponent);
		}
	private:
		const double exponent;
};

// A soft potential with a pentagon shape
class SoftPentagon : public PotentialType {
	public:
		SoftPentagon() {
			description = "softpentagon";
		}
		inline double operator()(double x, double y) const {
			const double r2 = x*x+y*y;
			const double t = atan2(y,x);
			return 0.5*r2*(1 + 0.5*sin(5*t));
		}
};

// A Hénon-Heiles type potential, modified so that all trajectories are bounded and generalized
// as in http://mathworld.wolfram.com/Henon-HeilesEquation.html
class HenonHeiles : public PotentialType {
	public:
		HenonHeiles(double _a=205.0/42.0, double _b=-13.0/3.0) : a(_a), b(_b) {
			std::stringstream ss;
			ss << "henonheiles(" << a << "," << b << ")";
			description = ss.str();
		}
		inline double operator()(double x, double y) const {
			const double r2 = x*x+y*y;
			const double r3 = pow(r2, 3.0/2.0);
			const double r4 = r2*r2;
			const double t = atan2(y,x);
			return r4 + a*r2 + b*r3*cos(3*t);
		}
	private:
		const double a;
		const double b;
};

// A Gaussian blob with a prescribed amplitude and width, centered at (x0,y0)
class Gaussian : public PotentialType {
	public:
		Gaussian(double _amplitude=1, double _width=1, double _x0=0, double _y0=0) :
				amplitude(_amplitude), width(_width), x0(_x0), y0(_y0) {
			if (amplitude < 0)
				throw InvalidPotentialType("gaussian potential with negative amplitude");
			if (width == 0)
				throw InvalidPotentialType("gaussian potential with zero width");
			std::stringstream ss;
			ss << "gaussian(" << amplitude << "," << width << "," << x0 << "," << y0 << ")";
			description = ss.str();
		}
		inline double operator()(double x, double y) const {
			const double xp = x-x0;
			const double yp = y-y0;
			const double exponent = -(xp*xp + yp*yp)/(2*width*width);
			return amplitude*exp(exponent);
		}
	private:
		const double amplitude;
		const double width;
		const double x0;
		const double y0;
};

#endif // _POTENTIALTYPES_HPP_
