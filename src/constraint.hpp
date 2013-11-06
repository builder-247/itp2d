/* Copyright 2013 Perttu Luukko

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
 * A simple class representing a geometric constraint (a yes/no answer for each
 * point in space), for example in order to constraint noise added to a
 * potential to a certain area.
 */

#ifndef _CONSTRAINT_HPP_
#define _CONSTRAINT_HPP_

#include <cmath>
#include <vector>
#include <tr1/tuple>

#include "itp2d_common.hpp"
#include "parser.hpp"

/* Interface class representing a geometric constraint */
class Constraint {
	public:
		virtual ~Constraint() {};
		// The main constraint function: return true for allowed points and
		// false for disallowed points
		virtual bool check(double x, double y) const = 0;
		virtual std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

// A parser function for returning a Noise instance from a user provided
// desciption string

Constraint const* parse_constraint_description(std::string const& str);
Constraint const* parse_constraint_description(name_parameters_pair const& pair);

// The trivial case of no constraint at all

class NoConstraint : public Constraint {
	public:
		NoConstraint() { init(); }
		NoConstraint(std::vector<double> params);
		bool check(__attribute__((unused)) double x, __attribute__((unused)) double y) const { return true; }
	private:
		void init() { description = "none"; }
};


// A special constraint representing a negation of another constraint

class InverseConstraint : public Constraint {
	public:
		InverseConstraint(Constraint const& base);
		InverseConstraint(Constraint const* base); // If constructed with a pointer, this pointer will be
												   // free'd in the destructor
		~InverseConstraint();
		inline bool check(double x, double y) const { return !base_constraint.check(x, y); }
	private:
		Constraint const& base_constraint;
		Constraint const* base_constraint_ptr;
		const bool owns_base_constraint;
		void init();
};

// Maximum radial distance

class MaximumRadialDistanceConstraint : public Constraint {
	public:
		MaximumRadialDistanceConstraint(double _r) : r(_r) { init(); }
		MaximumRadialDistanceConstraint(std::vector<double> params);
		bool check(__attribute__((unused)) double x, __attribute__((unused)) double y) const {
			return (hypot(x,y) <= r);
		}
	private:
		double r;
		void init();
};

// A ring-like (or annulus-like) constraint with a minimum radius and a width

class RingConstraint : public Constraint {
	public:
		RingConstraint(double _minr, double _width) : minr(_minr), width(_width) { init(); }
		RingConstraint(std::vector<double> params);
		bool check(__attribute__((unused)) double x, __attribute__((unused)) double y) const {
			const double r = hypot(x,y);
			return (r >= minr and r <= minr+width);
		}
	private:
		double minr;
		double width;
		void init();
};

#endif // _CONSTRAINT_HPP_
