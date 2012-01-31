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

#include "operatorproduct.hpp"

std::ostream& OperatorProduct::print(std::ostream& stream) const {
	for (const_opiter op = components.begin(); op != components.end(); op++) {
		stream << (**op);
	}
	return stream;
}

OperatorProduct::OperatorProduct() : components() {}

OperatorProduct::OperatorProduct(Operator const& oper) : components(1,&oper) {}

OperatorProduct::OperatorProduct(OperatorProduct const& operprod)
	: components(operprod.components.begin(), operprod.components.end()) {}

// The size required by an OperatorProduct is simply the maximum of the sizes
// required by its components
size_t OperatorProduct::required_workspace() const {
	size_t size = 0;
	for (const_opiter op = components.begin(); op != components.end(); op++) {
		const size_t req = (**op).required_workspace();
		if (req > size)
			size = req;
	}
	return size;
}

// OperatorProduct will simply act on the state with all of its components in order
void OperatorProduct::operate(State& state, StateArray& workspace) const {
	for (const_ropiter op = components.rbegin(); op != components.rend(); op++) {
		(**op)(state, workspace);
	}
}
