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

#ifndef _CONVERGENCE_HPP_
#define _CONVERGENCE_HPP_

/* Classes representing different tests for convergence */

#include "itp2d_common.hpp"

// Forward declaration of class ITPSystem; actual implementation will be included in the cpp file
class ITPSystem;

// Main interface class

class ConvergenceTest {
	public:
		// Main function for testing whether state 'n' is converged according to the test.
		// Since convergence testing can be very general, we pass the function the whole ITPSystem data (read-only).
		// This could be done with better data isolation, but let this be an exercise for the reader.
		virtual ~ConvergenceTest() {}
		virtual bool test(ITPSystem const& sys, size_t n) const = 0;	// returns true if state 'n' is converged in system sys.
		inline std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

// Convergence tests

// A simple no-op test that always fails, i.e., states are propagated until the maximum iteration count is reached.
class NoConvergenceTest : public ConvergenceTest {
	public:
		NoConvergenceTest() {
			description = "none";
		}
		bool test(__attribute__((unused)) ITPSystem const& sys, __attribute__((unused)) size_t n) const { return false; }
};

// Consider a state completely converged when it reaches convergence w.r.t
// current ltimestep size with one iteration.
class OneStepConvergenceTest : public ConvergenceTest {
	public:
		OneStepConvergenceTest() {
			description = "one-step convergence";
		}
		bool test(ITPSystem const& sys, size_t n) const;
};

// Consider a state converged if the relative (or absolute) change in energy
// between successive steps is less than a specified limit.
class RelativeEnergyChangeTest : public ConvergenceTest {
	public:
		RelativeEnergyChangeTest(double _limit) : limit(_limit) {
			std::stringstream ss;
			ss << "relative energy change < " << limit;
			description = ss.str();
		}
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		const double limit;
};

class AbsoluteEnergyChangeTest : public ConvergenceTest {
	public:
		AbsoluteEnergyChangeTest(double _limit) : limit(_limit) {
			std::stringstream ss;
			ss << "absolute energy change < " << limit;
			description = ss.str();
		}
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		const double limit;
};

// Calculate the variance of the state w.r.t the Hamiltonian, i.e., calculate
// <p|H²|p> - <p|H|p>², where H is the Hamiltonian and |p> is the state. If the
// square root of the variance -- the standard deviation -- is below
// a certain limit, consider the state converged. Also consider the state
// converged if the change in this standard deviation is small enough between
// successive steps, since this signals that the standard deviation is as small
// as it will get.
class EnergyDeviationChangeTest : public ConvergenceTest {
	public:
		EnergyDeviationChangeTest(double _relative_deviation_limit, double _difference_limit = 0) :
				relative_deviation_limit(_relative_deviation_limit),
				difference_limit(_difference_limit) {
			std::stringstream ss;
			ss << "relative energy deviation < " << relative_deviation_limit
				<< " or relative energy deviation change < " <<
				difference_limit;
			description = ss.str();
		}
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		const double relative_deviation_limit;
		const double difference_limit;
};

#endif // _CONVERGENCE_HPP_
