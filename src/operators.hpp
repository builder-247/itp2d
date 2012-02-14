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

#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "itp2d_common.hpp"
#include <ostream>
#include <string>

#include "state.hpp"
#include "statearray.hpp"

class Operator {
	public:
		virtual inline ~Operator() {};
		virtual void operate(State& state, __attribute__((unused)) StateArray& workspace) const = 0;
		virtual size_t required_workspace() const = 0;	// The required size (in number of copies of State) for operator()
		inline void operator()(State& state, StateArray& workspace) const;
		virtual std::ostream& print(std::ostream& out) const = 0;
		// The matrix element <p|O|s>, where |p> is the left state, |s> is the right state and O is this operator.
		// We can also have an exponent e, in which case we calculate <p|O^e|s>.
		comp matrixelement(State const& left, State const& right, StateArray& workspace, int exponent=1) const;
		// The expected value and standard deviation of the operator when acting on a state.
		std::pair<comp,comp> mean_and_standard_deviation(State const& state, StateArray& workspace) const;
		inline comp standard_deviation(State const& state, StateArray& workspace) const;
};

inline void Operator::operator()(State& state, StateArray& workspace) const {
			assert(state.datalayout == workspace.datalayout);
			assert(workspace.size() >= required_workspace());
			operate(state, workspace);
}

inline comp Operator::standard_deviation(State const& state, StateArray& workspace) const {
	const std::pair<comp,comp> masd = mean_and_standard_deviation(state, workspace);
	return masd.second;
}

std::ostream& operator<<(std::ostream& stream, __attribute__((unused)) const Operator& op);

// Simple interface class for evolution operators. These have the common
// property that they are dependent on the imaginary time step value, which can
// change.

class EvolutionOperator : public virtual Operator {
	public:
		inline virtual ~EvolutionOperator() {};
		virtual void set_time_step(double e) = 0;
};

#endif // _OPERATORS_HPP_
