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
 * A class representing a complete ITP simulation hiding as much internal
 * details as possible under the hood.
 */

#ifndef _ITPSYSTEM_HPP_
#define _ITPSYSTEM_HPP_

#include <string>
#include <ostream>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <ctime>
#include <csignal>
#include <tr1/tuple>

#include <omp.h>

#include "itp2d_common.hpp"
#include "exceptions.hpp"
#include "datafile.hpp"
#include "state.hpp"
#include "stateset.hpp"
#include "statearray.hpp"
#include "operators.hpp"
#include "potential.hpp"
#include "kinetic.hpp"
#include "multiproductsplit.hpp"
#include "rng.hpp"
#include "timer.hpp"
#include "parameters.hpp"
#include "potentialtypes.hpp"
#include "convergence.hpp"

class ITPSystem {
	public:
		typedef std::tr1::tuple<double,double,size_t> Esn_tuple;
		// Constructors & destructors
		// All parameters for ITPSystem are provided by the Parameters class
		ITPSystem(Parameters const& params,
				volatile sig_atomic_t* abort_flagptr = NULL,
				volatile sig_atomic_t* save_flagptr = NULL,
				std::ostream& out = std::cout,
				std::ostream& err = std::cerr);
		~ITPSystem();
		// Status checks
		inline size_t how_many_timestep_converged() { return states.get_num_timestep_converged(); }
		inline size_t how_many_finally_converged() { return states.get_num_finally_converged(); }
		size_t lowest_not_finally_converged();
		inline bool get_error_flag() const { return error_flag; }
		inline bool is_finished() const { return finished; }
		inline int get_total_step_counter() const { return total_step_counter; }
		inline int get_step_counter() const { return step_counter; }
		// Getters
		inline std::vector<std::vector<double> > const& get_energies() const { return energies; }
		inline std::vector<std::vector<double> > const& get_standard_deviations() const { return standard_deviations; }
		inline double get_sorted_energy(size_t n) const { return std::tr1::get<0>(Esn_tuples[n]); }
		inline size_t get_sorted_index(size_t n) const { return std::tr1::get<2>(Esn_tuples[n]); }
		inline double get_total_time() { return total_timer.get_time(); }
		inline double get_prop_time() { return prop_timer.get_time(); }
		inline double get_ortho_time() { return states.get_ortho_time(); }
		inline double get_dot_time() { return states.get_dot_time(); }
		inline double get_eigensolve_time() { return states.get_eigensolve_time(); }
		inline double get_lincomb_time() { return states.get_lincomb_time(); }
		inline double get_io_time() { return io_timer.get_time(); }
		inline double get_convtest_time() { return convtest_timer.get_time(); }
		inline StateSet const& get_states() const { return states; }
		inline State const& get_state(size_t n) const { return states[n]; }
		inline Potential const& get_potential() const { return *pot; }
		inline OperatorSum const& get_hamiltonian() const { return H; }
		inline double get_eps() const { return eps; }
		// Main operation
		void step();	// A single iteration of ITP
		void check_timestep_convergence();
		void check_final_convergence();
		void calculate_energies();
		void save_states(bool sort = true); // If sort is true, states are sorted according to energy
		void save_energies();
		void save_energy_history();
		void print_energies();
		void finish();
		const Parameters params;
		const DataLayout datalayout;
		const Transformer transformer;
		const BoundaryType boundary_type;
	private:
		void print_initial_message();
		void print_final_message();
		void propagate();
		void orthonormalize();
		void change_time_step();
		inline void check_save_flag();
		inline bool verb(int level) const { return (params.get_verbosity() >= level)? true : false; }
		inline void update_timestring();
		// Pointers to flags set by external signals
		volatile sig_atomic_t* abort_flagptr;
		volatile sig_atomic_t* save_flagptr;
		// Output streams
		std::ostream& out;
		std::ostream& err;
		// Status flags
		bool finished;
		bool error_flag;
		bool all_needed_states_timestep_converged;
		bool all_needed_states_finally_converged;
		bool exhausting_eps_values;
		// Timers, RNG etc helpers
		RNG rng;
		Timer total_timer, prop_timer, io_timer, convtest_timer;
		time_t rawtime;
		struct tm* timeinfo;
		char timestring[24];
		// Main members
		PotentialType const* pot_type;
		Noise const* noise;
		ImpurityType const* impurity_type;
		ImpurityDistribution const* impurity_distribution;
		Constraint const* impurity_constraint;
		Potential const* pot;
		const Kinetic kin;
		OperatorSum H;
		MultiProductSplit* T;
		Datafile* datafile;
		StateSet states;
		StateArray** workslices;
		std::vector<std::vector<double> > energies;				// A vector of energy values for each iteration
		std::vector<std::vector<double> > standard_deviations;	// ... and the same thing for the standard deviations of energy
		std::vector<Esn_tuple> Esn_tuples;	// A vector of tuples (E,s,n), where E is the energy of a state,
											// s is the standard deviation, and n is the index where the state is stored in the StateSet.
											// This will be updated whenever new energy values are calculated.
		// Running counters etc.
		int total_step_counter;
		int step_counter;
		double eps;
		std::list<double> eps_values;
};

inline void ITPSystem::update_timestring() {
	time (&rawtime);
	timeinfo = gmtime(&rawtime);
	__attribute__((unused)) const size_t retval = strftime(timestring, 24, "%Y-%m-%d %H:%M:%SZ", timeinfo);
	assert(retval != 0);
}

#endif // _ITPSYSTEM_HPP_
