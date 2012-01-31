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

#include "secondordersplit.hpp"

SecondOrderSplit::SecondOrderSplit(Potential const& original_potential, double time_step, double B, Transformer const& tr, BoundaryType bt, double prefactor, int exponent) {
	assert(tr.datalayout == original_potential.datalayout);
	if (not original_potential.is_null()) {
		kinetic_part = new ExpKinetic(time_step, B, tr, bt, -1.0);
		potential_part = new ExpPotential(original_potential, time_step, -0.5, 1.0);
		if (prefactor != 1.0)
			potential_with_prefactor = new ExpPotential(original_potential, time_step, -0.5, prefactor);
		else
			potential_with_prefactor = potential_part;
		if (exponent >= 2)
			potential_part_square = new ExpPotential(original_potential, time_step, -1.0, 1.0);
		else
			potential_part_square = NULL;
		(*this) *= (*potential_with_prefactor);
		(*this) *= (*kinetic_part);
		for (int t=2; t<=exponent; t++) {
			(*this) *= (*potential_part_square);
			(*this) *= (*kinetic_part);
		}
		(*this) *= (*potential_part);
	}
	else { // No splitting needed if potential is zero
		kinetic_part = new ExpKinetic(time_step, B, tr, bt, -1.0*exponent, prefactor);
		(*this) *= (*kinetic_part);
		potential_part = NULL;
		potential_part_square = NULL;
		potential_with_prefactor = NULL;
	}
}

SecondOrderSplit::~SecondOrderSplit() {
	delete kinetic_part;
	delete potential_part_square;
	// Avoid double free
	if (potential_with_prefactor != potential_part) {
		delete potential_with_prefactor;
	}
	delete potential_part;
}

void SecondOrderSplit::set_time_step(double time_step) {
	kinetic_part->set_time_step(time_step);
	if (potential_part != NULL)
		potential_part->set_time_step(time_step);
	if (potential_part_square != NULL)
		potential_part_square->set_time_step(time_step);
	if (potential_with_prefactor != potential_part and potential_with_prefactor != NULL)
		potential_with_prefactor->set_time_step(time_step);
}
