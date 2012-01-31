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

#ifndef _OPERATORPRODUCT_HPP_
#define _OPERATORPRODUCT_HPP_

/*
 *	A product of operators is simply an operator which applies each operator in
 *	the chain and feeds the result into the next operator.
 */

#include <list>
#include "operators.hpp"

class OperatorProduct : public virtual Operator {
	public:
		// helpful typedefs
		typedef std::list<Operator const*> oplist;
		typedef oplist::iterator opiter;
		typedef oplist::const_iterator const_opiter;
		typedef oplist::reverse_iterator ropiter;
		typedef oplist::const_reverse_iterator const_ropiter;
		// the empty product
		OperatorProduct();
		// a product of a single operator
		OperatorProduct(Operator const& oper);
		// copy constructor
		OperatorProduct(OperatorProduct const& operprod);
		void operate(State& state, StateArray& workspace) const;
		size_t required_workspace() const;
		std::ostream& print(std::ostream& out) const;
		// Arithmetic
		inline OperatorProduct& operator*=(Operator const& oper) { components.push_back(&oper); return *this; }
		inline OperatorProduct& operator*=(OperatorProduct const& operprod);
	protected:
		oplist components;
};

inline OperatorProduct& OperatorProduct::operator*=(OperatorProduct const& operprod) {
	for (const_ropiter op = operprod.components.rend(); op != operprod.components.rbegin(); op++) {
		components.push_back(*op);
	}
	return *this;
}

// Free functions allowing Operators to be multiplied simply with *

inline OperatorProduct operator*(Operator const& lhand, Operator const& rhand) {
	return OperatorProduct(lhand) *= rhand;
}

inline OperatorProduct operator*(OperatorProduct const& lhand, Operator const& rhand) {
	return OperatorProduct(lhand) *= rhand;
}

inline OperatorProduct operator*(Operator const& lhand, OperatorProduct const& rhand) {
	return OperatorProduct(lhand) *= rhand;
}

#endif // _OPERATORPRODUCT_HPP_
