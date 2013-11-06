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

#include "constraint.hpp"

// Constraint parser function

Constraint const* parse_constraint_description(std::string const& str) {
	// Parse with the generic parse_parameter_string function and extract the
	// name and parameters.
	name_parameters_pair p;
	try {
		p = parse_parameter_string(str);
	}
	catch (ParseError& e) {
		std::cerr << e.what() << std::endl;
		throw UnknownConstraintType(str);
	}
	return parse_constraint_description(p);
}

Constraint const* parse_constraint_description(name_parameters_pair const& p) {
	std::string const& name = p.first;
	std::vector<double> const& params = p.second;
	// Simply delegate to the individual constructors based on name
	if (name[0] == '!') {
		std::string base = name.substr(1,name.size()-1);
		return new InverseConstraint(parse_constraint_description(make_pair(base, params)));
	}
	else if (name == "no" or name == "none" or name == "zero")
		return new NoConstraint(params);
	else if (name == "maxr" or name == "maxradius")
		return new MaximumRadialDistanceConstraint(params);
	else if (name == "ring" or name == "annulus")
		return new RingConstraint(params);
	else
		throw UnknownConstraintType(name);
}

// NoConstraint

NoConstraint::NoConstraint(std::vector<double> params) {
	if (not params.empty())
		throw InvalidConstraintType("Constraint type NoConstraint does not take parameters");
	init();
}

// InverseConstraint

InverseConstraint::InverseConstraint(Constraint const& base) :
		base_constraint(base),
		base_constraint_ptr(&base),
		owns_base_constraint(false) {
	init();
}

InverseConstraint::InverseConstraint(Constraint const* base) :
		base_constraint(*base),
		base_constraint_ptr(base),
		owns_base_constraint(true) {
	init();
}

InverseConstraint::~InverseConstraint() {
	if (owns_base_constraint)
		delete base_constraint_ptr;
}

void InverseConstraint::init() {
	std::stringstream ss;
	ss << "inverse of " << base_constraint.get_description();
	description = ss.str();
}

// MaximumRadialDistanceConstraint

MaximumRadialDistanceConstraint::MaximumRadialDistanceConstraint(std::vector<double> params) {
	if (params.size() == 1) {
		r = params[0];
	}
	else
		throw InvalidConstraintType("Constraint type MaximumRadialDistanceConstraint takes exactly 1 parameter");
	init();
}

void MaximumRadialDistanceConstraint::init() {
	if (r < 0)
		throw InvalidConstraintType("radial constraint with negative radius");
	std::stringstream ss;
	ss << "maximum radial distance = " << r;
	description = ss.str();
}

// RingConstraint

RingConstraint::RingConstraint(std::vector<double> params) {
	if (params.size() == 2) {
		minr = params[0];
		width = params[1];
	}
	else
		throw InvalidConstraintType("Constraint type RingConstraint takes exactly 2 parameters");
	init();
}

void RingConstraint::init() {
	if (minr < 0)
		throw InvalidConstraintType("ring constraint with negative inner radius");
	if (width < 0)
		throw InvalidConstraintType("ring constraint with negative width");
	std::stringstream ss;
	ss << "ring with inner radius = " << minr << " and width = " << width;
	description = ss.str();
}
