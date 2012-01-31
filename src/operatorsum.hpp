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

#ifndef _OPERATORSUM_HPP_
#define _OPERATORSUM_HPP_

/* 
 * A sum of operators is an operator which applies each operator in the chain to the same state and adds the results together.
 */

#include <list>
#include "operators.hpp"

class OperatorSum : public virtual Operator {
	public:
		// helpful typedefs first
		typedef std::list<Operator const*> oplist;
		typedef oplist::iterator opiter;
		typedef oplist::const_iterator const_opiter;
		typedef oplist::reverse_iterator ropiter;
		typedef oplist::const_reverse_iterator const_ropiter;
		OperatorSum();
		OperatorSum(Operator const& oper);
		OperatorSum(OperatorSum const& opersum);
		void operate(State& state, StateArray& workspace) const;
		size_t required_workspace() const;
		std::ostream& print(std::ostream& out) const;
		// Arithmetic
		inline OperatorSum& operator+=(Operator const& oper) { components.push_back(&oper); return *this; }
		inline OperatorSum& operator+=(OperatorSum const& operprod);
	protected:
		oplist components;
};

inline OperatorSum& OperatorSum::operator+=(OperatorSum const& opersum) {
	for (const_opiter op = opersum.components.end(); op != opersum.components.begin(); op++) {
		components.push_back(*op);
	}
	return *this;
}

// Free functions which allow operators to be summed simply with +

inline OperatorSum operator+(Operator const& lhand, Operator const& rhand) {
	return OperatorSum(lhand) += rhand;
}

inline OperatorSum operator+(OperatorSum const& lhand, Operator const& rhand) {
	return OperatorSum(lhand) += rhand;
}

inline OperatorSum operator+(Operator const& lhand, OperatorSum const& rhand) {
	return OperatorSum(lhand) += rhand;
}

#endif // _OPERATORSUM_HPP_
