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
#include "parser.hpp"

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

// Parser function for returning a ConvergenceTest* from a user-provided string

ConvergenceTest const* parse_convergence_description(std::string const& str);

// Convergence tests

// A simple no-op test that always fails, i.e., states are propagated until the maximum iteration count is reached.
class NoConvergenceTest : public ConvergenceTest {
	public:
		NoConvergenceTest() { init(); }
		NoConvergenceTest(std::vector<double> const& params);
		bool test(__attribute__((unused)) ITPSystem const& sys, __attribute__((unused)) size_t n) const { return false; }
	private:
		void init() { description = "none"; }
};

// Consider a state completely converged when it reaches convergence w.r.t
// current timestep size with one iteration.
class OneStepConvergenceTest : public ConvergenceTest {
	public:
		OneStepConvergenceTest() { init(); }
		OneStepConvergenceTest(std::vector<double> const& params);
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		void init() { description = "one-step convergence"; }
};

// Consider a state converged if the relative (or absolute) change in energy
// between successive steps is less than a specified limit.
class RelativeEnergyChangeTest : public ConvergenceTest {
	public:
		RelativeEnergyChangeTest(double _limit) : limit(_limit) { init(); }
		RelativeEnergyChangeTest(std::vector<double> const& params);
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		double limit;
		void init();
};

class AbsoluteEnergyChangeTest : public ConvergenceTest {
	public:
		AbsoluteEnergyChangeTest(double _limit) : limit(_limit) { init(); }
		AbsoluteEnergyChangeTest(std::vector<double> const& params);
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		double limit;
		void init();
};

// Calculate the variance of the state w.r.t the Hamiltonian, i.e., calculate
// <p|H²|p> - <p|H|p>², where H is the Hamiltonian and |p> is the state. If the
// square root of the variance -- the standard deviation -- is below
// a certain limit, consider the state converged. Also consider the state
// converged if the change in this standard deviation is small enough between
// successive steps, since this signals that the standard deviation is as small
// as it will get.
// This is also the standard residual value used in eigenvalue solving.
class AbsoluteEnergyDeviationTest : public ConvergenceTest {
	public:
		AbsoluteEnergyDeviationTest(double _absolute_deviation_limit, double _difference_limit = 0) :
				absolute_deviation_limit(_absolute_deviation_limit),
				difference_limit(_difference_limit) {
			init();
		}
		AbsoluteEnergyDeviationTest(std::vector<double> const& params);
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		double absolute_deviation_limit;
		double difference_limit;
		void init();
};

// The same as the previous case, but use a relative standard deviation.
class RelativeEnergyDeviationTest : public ConvergenceTest {
	public:
		RelativeEnergyDeviationTest(double _relative_deviation_limit, double _difference_limit = 0) :
				relative_deviation_limit(_relative_deviation_limit),
				difference_limit(_difference_limit) {
			init();
		}
		RelativeEnergyDeviationTest(std::vector<double> const& params);
		bool test(ITPSystem const& sys, size_t n) const;
	private:
		double relative_deviation_limit;
		double difference_limit;
		void init();
};

#endif // _CONVERGENCE_HPP_
