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
		virtual inline ~PotentialType() {}
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
		ZeroPotential() { init(); }
		ZeroPotential(std::vector<double> params);
		inline double operator()(__attribute__((unused)) double x, __attribute__((unused)) double y) const { return 0; }
	private:
		void init() { description = "zero"; }
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
class HarmonicOscillator : public PotentialType {
	public:
		static const double default_frequency;
		static const double default_x0;
		static const double default_y0;
		HarmonicOscillator(double omega=default_frequency, double orig_x=default_x0, double orig_y=default_y0);
		HarmonicOscillator(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double px = x-x0;
			const double py = y-y0;
			return 0.5*w*(px*px+py*py);
		}
	private:
		double w;
		double x0;
		double y0;
		void init();
};

// The non-degenerate (elliptic) harmonic oscillator
class EllipticOscillator : public PotentialType {
	public:
		static const double default_frequency_x;
		static const double default_frequency_y;
		static const double default_x0;
		static const double default_y0;
		EllipticOscillator(double omega_x=default_frequency_x, double omega_y=default_frequency_y);
		EllipticOscillator(std::vector<double> params);
		inline double operator()(double x, double y) const {
			return 0.5*(wx*x*x + wy*y*y);
		}
	private:
		double wx;
		double wy;
		void init();
};

// A pretty hard square of length pi
class PrettyHardSquare : public PotentialType {
	public:
		static const double default_exponent;
		PrettyHardSquare(double e=default_exponent) : exponent(e) { init(); }
		PrettyHardSquare(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double ax = fabs(2*x/pi);
			const double ay = fabs(2*y/pi);
			const double r = (ax > ay)? ax : ay;
			return pow(r,exponent);
		}
	private:
		double exponent;
		void init();
};

// A soft potential with a pentagon shape
class SoftPentagon : public PotentialType {
	public:
		SoftPentagon() { init(); }
		SoftPentagon(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double r2 = x*x+y*y;
			const double t = atan2(y,x);
			return 0.5*r2*(1 + 0.5*sin(5*t));
		}
	private:
		void init() { description = "softpentagon"; }
};

// A Hénon-Heiles type potential, modified so that all trajectories are bounded and generalized
// as in http://mathworld.wolfram.com/Henon-HeilesEquation.html
class HenonHeiles : public PotentialType {
	public:
		static const double default_a;
		static const double default_b;
		HenonHeiles(double _a=default_a, double _b=default_b) : a(_a), b(_b) { init(); }
		HenonHeiles(std::vector<double> params);
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
		void init();
};

// A GaussianPotential blob with a prescribed amplitude and width, centered at (x0,y0)
class GaussianPotential : public PotentialType {
	public:
		static const double default_amplitude;
		static const double default_width;
		static const double default_x0;
		static const double default_y0;
		GaussianPotential(double _amplitude=default_amplitude, double _width=default_width,
				double _x0=default_x0, double _y0=default_y0) :
				amplitude(_amplitude), width(_width), x0(_x0), y0(_y0) {
			init();
		}
		GaussianPotential(std::vector<double> params);
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
		void init();
	};

// The quartic oscillator, frequently studied in quantum chaos.
// The potential is P(x,y) = (x^2 * y^2)/2 + b(x^4 + y^4)/4,
// but it is rotated by pi/4 to make maximal use of the rectangular
// grid.
// For references please see e.g.
// 	* Eckhardt et al, Phys Rev A 39 3776 (1989)
//	* de Polavieja et al, Phys Rev Lett 73 1613 (1994)
class QuarticPotential : public PotentialType {
	public:
		static const double default_b;
		QuarticPotential(double _b = default_b) : b(_b) { init(); }
		QuarticPotential(std::vector<double> params);
		inline double operator()(double xp, double yp) const {
			// rotate first
			const double x = (xp - yp)/sqrt(2);
			const double y = (xp + yp)/sqrt(2);
			// then use the simple formula
			const double x2 = x*x;
			const double y2 = y*y;
			const double x4 = x2*x2;
			const double y4 = y2*y2;
			return (x2*y2)/2 + b*(x4+y4)/4;
		}
	private:
		double b;
		void init();
};

//
class SquareOscillator : public PotentialType {
	public:
		static const double default_alpha;
		SquareOscillator(double alpha_=default_alpha) : alpha(alpha_) { init(); }
		SquareOscillator(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double ax = fabs(x);
			const double ay = fabs(y);
			return 0.5*(pow(ax,alpha) + pow(ay,alpha));
		}
	private:
		double alpha;
		void init();
};

// A generalization of the harmonic oscillator
class PowerOscillator : public PotentialType {
	public:
		static const double default_exponent;
		static const double default_w;
		PowerOscillator(double exponent=default_exponent, double omega=default_w);
		PowerOscillator(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double r = hypot(x, y);
			return 0.5*w*pow(r,a);
		}
	private:
		double a;
		double w;
		void init();
};

// A flexible class for ring-like potentials, with optional asymmetry in the
// form of a Gaussian spike
class RingPotential : public PotentialType {
	public:
		static const double default_radius;
		static const double default_width;
		static const double default_exponent;
		static const double default_asymm_amplitude;
		static const double default_asymm_width;
		RingPotential(double radius=default_radius, double width=default_width, double exponent=default_exponent,
				double asymm_amplitude=default_asymm_amplitude, double asymm_width=default_asymm_width);
		RingPotential(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double rp = hypot(x, y);
			const double z = fabs(rp-r)/w;
			const double xp = x-r;
			const double G_exponent = -(xp*xp + y*y)/(2*asymm_w*asymm_w);
			const double G = asymm_A*exp(G_exponent);
			return 0.5*pow(z,e) + G;
		}
	private:
		double r;
		double w;
		double e;
		double asymm_A;
		double asymm_w;
		void init();
};

// A radial cosh-potential
class CoshPotential : public PotentialType {
	public:
		static const double default_amplitude;
		static const double default_length_scale;
		CoshPotential(double amplitude=default_amplitude, double length_scale=default_length_scale);
		CoshPotential(std::vector<double> params);
		inline double operator()(double x, double y) const {
			const double r = hypot(x, y);
			return A*(cosh(r/L)-1);
		}
	private:
		double A;
		double L;
		void init();
};

// Soft stadium potential as described by
// Tomsovic, S. and Heller, E. J., PRE 47, 282 (1993)
class SoftStadium : public PotentialType {
	public:
		static const double default_radius;
		static const double default_center_length;
		static const double default_height;
		static const double default_a;
		static const double default_b;
		SoftStadium(double radius=default_radius, double center_length=default_center_length,
				double height=default_height, double a=default_a, double b=default_b);
		SoftStadium(std::vector<double> params);
		inline double operator()(double x, double y) const {
			double q;
			if (fabs(x) <= halfL) { // in the central rectangle
				q = y/R;
			}
			else if (x < -halfL) { //  in the left end cap
				q = hypot(x+halfL, y)/R;
			}
			else { // in the right end cap
				q = hypot(x-halfL, y)/R;
			}
			return V/(1+a*exp(b*(1-q*q)));
		}
	private:
		double R;
		double halfL;
		double V;
		double a;
		double b;
		void init();
};

#endif // _POTENTIALTYPES_HPP_
