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

#include "exppotential.hpp"

std::ostream& ExpPotential::print(std::ostream& stream) const {
	if (prefactor != 1.0)
		stream << prefactor << "·";
	stream << "exp(" << time_step*coefficient << "·"  << original_name << ")";
	return stream;
}

ExpPotential::ExpPotential(Potential const& pot, double e, double c, double p) :
		datalayout(pot.datalayout),
		original_potential(pot),
		values(NULL),
		time_step(e),
		coefficient(c),
		prefactor(p),
		original_name(pot.get_name()),
		is_trivial(pot.is_null()) {
	if (not is_trivial) {	// Don't bother calculating an array of ones if the potential is always zero.
		values = new double[datalayout.N];
		recalc_potential();
	}
}

ExpPotential::~ExpPotential() {
	if (not is_trivial)
		delete[] values;
}

void ExpPotential::set_constants(double e, double c, double p) {
	time_step = e;
	coefficient = c;
	prefactor = p;
	recalc_potential();
}

void ExpPotential::recalc_potential() {
	if (is_trivial)
		return;
	for (size_t x=0; x<datalayout.sizex; x++) {
		for (size_t y=0; y<datalayout.sizey; y++) {
			datalayout.value(values, x, y) =  prefactor*exp(time_step*coefficient * original_potential.get_value(x,y));
		}
	}
}
