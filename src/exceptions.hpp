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

#ifndef _EXCEPTIONS_HPP_
#define _EXCEPTIONS_HPP_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>

#include "itp2d_common.hpp"

// A helper function for pretty-printing a matrix for some seriously deep error reporting.
void dump_matrix(std::ostream& s, comp const* M, size_t N);

// Exception classes so that we have something to throw when things go wrong.
// These are all simply versions of std::runtime_error with different strings
// to show to the user.

class GeneralError : public std::runtime_error {
	public:
		GeneralError(std::string whatiswrong) : std::runtime_error("") {
			std::stringstream ss;
	        ss << "Error: " << whatiswrong;
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class NotImplemented : public std::runtime_error {
	public:
		NotImplemented(std::string missingfeature) : std::runtime_error("") {
			std::stringstream ss;
	        ss << missingfeature << " has not been implemented yet.";
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class NonPositiveEigenvalue : public std::runtime_error {
	public:
		NonPositiveEigenvalue(size_t n, double d, __attribute__((unused)) comp const* M, __attribute__((unused)) size_t N) : std::runtime_error("") {
			std::stringstream ss;
	        ss << "Got a non-positive eigenvalue " << d << " in orthonormalization at position " << n << "." << std::endl
				<< "This means that the states that were being orthonormalized turned out to be linearly dependent," << std::endl
				<< "or at least close enough. This can happen if the states are propagated \"too much\" in one go." << std::endl
				<< "Try with a smaller value of the time step.";
			#ifndef NDEBUG
			if (M != NULL and N <= 8) {
				ss << std::endl << "Raw matrix:" << std::endl;
				dump_matrix(ss, M, N);
			}
			#endif
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class NonNormalEigenvalue : public std::runtime_error {
	public:
		NonNormalEigenvalue(size_t n, double d, __attribute__((unused)) comp const* M, __attribute__((unused)) size_t N) : std::runtime_error("") {
			std::stringstream ss;
	        ss << "Got a non-normal eigenvalue " << d << " in orthonormalization at position " << n << ".";
			#ifndef NDEBUG
			if (M != NULL and N <= 8) {
				ss << std::endl << "Raw matrix:" << std::endl;
				dump_matrix(ss, M, N);
			}
			#endif
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class NaNError : public std::runtime_error {
	public:
		NaNError(const char* place) : std::runtime_error("") {
			std::stringstream ss;
	        ss << "A NaN bounced up where it really shouldn't: " << place << ".";
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class NonHermitianOverlapMatrix : public std::runtime_error {
	public:
		NonHermitianOverlapMatrix(size_t i, size_t j, comp z, comp w) : std::runtime_error("") {
			std::stringstream ss;
	        ss << "Overlap matrix is not Hermitian for some reason." << std::endl;
			ss << "state[" << i << "]·state[" << j << "] = " << z << std::endl;
			ss << "state[" << j << "]·state[" << i << "] = " << w << std::endl;
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class InvalidConvergenceHistoryAccess : public std::runtime_error {
	public:
		InvalidConvergenceHistoryAccess(size_t n, size_t step, size_t history) : std::runtime_error("") {
			std::stringstream ss;
			ss << "StateSet::log_eigen_deviation called with invalid arguments (" << n
				<< ", " << step << ", " << history << ").";
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class UnknownPotentialType : public std::runtime_error {
	public:
		UnknownPotentialType(std::string str) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Potential type \"" << str << "\" not understood by parser.";
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class UnknownConvergenceType : public std::runtime_error {
	public:
		UnknownConvergenceType(std::string str) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Convergence test description \"" << str << "\" not understood by parser.";
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class UnknownNoiseType : public std::runtime_error {
	public:
		UnknownNoiseType(std::string str) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Noise description \"" << str << "\" not understood by parser.";
			static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class InvalidPotentialType : public std::runtime_error {
	public:
		InvalidPotentialType(std::string str) : std::runtime_error("") {
			std::stringstream ss;
			ss	<< "Invalid potential: " << str << std::endl
				<< "Remember that potentials should have no eigenstates with negative energy." << std::endl;
    	    static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class InvalidDataspaceSelection : public std::runtime_error {
	public:
		InvalidDataspaceSelection(size_t const* bbox) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Invalid dataspace selection with bounding box ("
				<< bbox[0] << "," << bbox[1] << "," << bbox[2] << "," << bbox[3] << ").";
			static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class EigensolverError : public std::runtime_error {
	public:
		EigensolverError(int errorcode) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Error in eigenvalue solving. ZHEEV reported error code " << errorcode << ".";
			static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

class InvalidNumberOfStates : public std::runtime_error {
	public:
		InvalidNumberOfStates(size_t N, size_t num_states, size_t ignore_lowest) : std::runtime_error("") {
			std::stringstream ss;
			ss << "Supplied total number of states " << N << " is not enough to converge "
				<< num_states << " states and ignore lowest " << ignore_lowest << " states.";
			static_cast<std::runtime_error&>(*this) = std::runtime_error(ss.str());
		}
};

#endif // _EXCEPTIONS_HPP_
