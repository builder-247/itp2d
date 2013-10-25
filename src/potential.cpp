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

#include "potential.hpp"

Potential::Potential(DataLayout const& dl, PotentialType const& ptype, std::string arg_name) :
		datalayout(dl), type(ptype), name(arg_name) {
	if (typeid(ptype) == typeid(ZeroPotential)) {
		isnull = true;
		values = NULL;
	}
	else {
		isnull = false;
		init_values();
	}
}

Potential::Potential(DataLayout const& dl, PotentialType const& ptype, RNG& rng, Noise const& noise, Constraint const& noise_constraint, std::string arg_name) :
		datalayout(dl), type(ptype), name(arg_name) {
	if (typeid(ptype) == typeid(ZeroPotential) and typeid(noise) == typeid(NoNoise)) {
		isnull = true;
		values = NULL;
	}
	else {
		isnull = false;
		init_values();
		noise.add_noise(datalayout, values, rng, noise_constraint);
	}
}

Potential::~Potential() {
	delete[] values;
}

std::ostream& Potential::print(std::ostream& stream) const {
	return stream << "V";
}

void Potential::init_values() {
	values = new double[datalayout.sizex*datalayout.sizey];
	for (size_t y=0; y<datalayout.sizey; y++) {
		double dy = datalayout.get_posy(y);
		for (size_t x=0; x<datalayout.sizex; x++) {
			double dx = datalayout.get_posx(x);
			datalayout.value(values, x, y) = type(dx,dy);
		}
	}
}
