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

#include "operatorsum.hpp"

std::ostream& OperatorSum::print(std::ostream& stream) const {
	const_opiter next_op = components.begin();
	++next_op;
	for (const_opiter op = components.begin(); op != components.end(); ++op) {
		stream << (**op);
		++next_op;
		if (next_op != components.end())
			stream << " + ";
	}
	return stream;
}

OperatorSum::OperatorSum() : components() {}

OperatorSum::OperatorSum(Operator const& oper) : components(1,&oper) {}

OperatorSum::OperatorSum(OperatorSum const& opersum)
	: components(opersum.components.begin(), opersum.components.end()) {}

size_t OperatorSum::required_workspace() const {
	// special cases for 0 or 1 operators
	if (components.empty())
		return 0;
	if (components.size() == 1)
		return components.front()->required_workspace();
	// then the general case
	size_t size = 0;
	for (const_opiter op = components.begin(); op != components.end(); ++op) {
		const size_t req = (*op)->required_workspace();
		if (req > size)
			size = req;
	}
	// we need two more:
	//	one for saving the original state
	//	one for the intermediate result
	return size + 2;
	// This size requirement could be lowered if wee peek into the operators
	// and apply all in-place operators sequentially
}

void OperatorSum::operate(State& state, StateArray& workspace) const {
	// special cases for 0 or 1 operators
	if (components.empty()) {
		state.zero();
		return;
	}
	if (components.size() == 1) {
		(*components.front())(state, workspace);
		return;
	}
	// then the general case
	State& orig = workspace[0];
	State& intermediate = workspace[1];
	orig = state;
	// The rest is working space for the components
	StateArray workslice = StateArray(workspace, 2);
	const_opiter op = components.begin();
	(**op)(state, workslice);	// Operate with first operator
	++op;
	while (op != components.end()) {
		intermediate = orig; // Load original state
		(**op)(intermediate, workslice); // Operate on it
		state += intermediate;  // Sum to state
		++op;
	}
}
