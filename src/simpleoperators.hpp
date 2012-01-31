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

#ifndef _SIMPLEOPERATORS_HPP_
#define _SIMPLEOPERATORS_HPP_

#include "operators.hpp"

// Simple operators mainly used to test the logic of OperatorSum and
// OperatorProduct in unit testing.

// A simple null operator class

class NullOperator : public virtual Operator {
	public:
		void operate(__attribute__((unused)) State& state, __attribute__((unused))StateArray& workspace) const {};
		inline size_t required_workspace() const { return 0; }
		std::ostream& print(std::ostream& out) const;
};

// A simple zero operator class

class ZeroOperator : public virtual Operator {
	public:
		void operate(State& state, __attribute__((unused))StateArray& workspace) const { state.zero(); }
		inline size_t required_workspace() const { return 0; }
		std::ostream& print(std::ostream& out) const;
};

// A simple constant operator class

class ConstantOperator : public virtual Operator {
	public:
		ConstantOperator(comp c);
		void operate(State& state, __attribute__((unused))StateArray& workspace) const { state *= multiplier; };
		inline size_t required_workspace() const { return 0; }
		inline comp const& get_multiplier() const { return multiplier; }
		std::ostream& print(std::ostream& out) const;
	private:
		comp multiplier;
};

#endif // _SIMPLEOPERATORS_HPP_
