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
 * Abstract classes representing operators acting on State objects.
 */

#include "operators.hpp"

std::ostream& operator<<(std::ostream& stream, __attribute__((unused)) const Operator& op) {
	return op.print(stream);
}

comp Operator::matrixelement(State const& left, State const& right, StateArray& workspace, int exponent) const {
	assert(workspace.size() >= 1+required_workspace());
	State& temp = workspace[0];
	temp = right;
	StateArray workslice = StateArray(workspace, 1);
	for (int i=0; i<exponent; i++) {
		(*this)(temp, workslice);
	}
	return left.dot(temp);
}

std::pair<comp,comp> Operator::mean_and_standard_deviation(State const& state, StateArray& workspace) const {
	assert(workspace.size() >= 1+required_workspace());
	StateArray workslice = StateArray(workspace, 1);
	State& temp = workspace[0];
	temp = state;
	(*this)(temp, workslice);
	const comp mean = state.dot(temp);
	(*this)(temp, workslice);
	const comp meansqr = state.dot(temp);
	return std::pair<comp,comp>(mean, sqrt(meansqr-mean*mean));
}
