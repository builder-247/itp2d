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
 * A simple wrapper class representing the parameters given to an ITP simulation.
 */

#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <list>
#include <string>
#include <iterator>
#include <cassert>

#include "itp2d_common.hpp"
#include "convergence.hpp"
#include "noise.hpp"
#include "constraint.hpp"
#include "transformer.hpp"

class Parameters {
	friend class CommandLineParser;
	public:
		// Enums
		enum SaveWhat { Nothing, OnlyEnergies, FinalStates, Everything };
		enum InitialStatePreset { UserSuppliedInitialState, CopyFromFile, Random };
		// Typedefs
		typedef double (*potfunc)(double, double);
		typedef comp (*initialstatefunc)(size_t, double, double);
		// Constructors & destructors
		Parameters();
		~Parameters();
		// Initialization
		inline void set_recover(bool val) { recover = val; }
		inline void set_random_seed(unsigned long int s) { rngseed = s; }
		void define_data_storage(std::string filename, SaveWhat save_what = FinalStates, bool clobber = false);
		void set_wisdom_file_name(std::string const& filename) { wisdom_file_name = filename; }
		void define_grid(size_t sizex, size_t sizey, double lenx, BoundaryType boundary = Periodic);
		void define_external_field(std::string const& ptype, double B=0);
		void define_initial_states(InitialStatePreset preset = Random);
		void define_initial_states(std::string description, initialstatefunc func);
		void set_num_states(size_t N, size_t wanted_to_converge, size_t ignore_lowest=default_ignore_lowest);
		inline void set_operator_splitting_halforder(int h) { halforder = h; }
		inline void add_eps_value(double e) { eps_values.push_back(e); }
		inline void set_time_step_divisor(double d) { eps_divisor = d; }
		inline void set_exhaust_eps(bool val) { exhaust_eps = val; }
		void set_bailout_limits(int max_steps, double min_time_step);
		inline void set_ortho_algorithm(OrthoAlgorithm alg) { ortho_alg = alg; }
		inline void set_timestep_convergence_test(ConvergenceTest* test) { timestep_convergence_test = test; }
		inline void set_timestep_convergence_test(std::string const& str);
		inline void set_final_convergence_test(ConvergenceTest* test) { final_convergence_test = test; }
		inline void set_final_convergence_test(std::string const& str);
		inline void set_noise_type(std::string const& str) { noise_type = str; }
		inline void set_noise_constraint_type(std::string const& str) { noise_constraint_type = str; }
		inline void set_fftw_flags(unsigned int fl) { fftw_flags = fl; }
		inline void set_verbosity(int val) { verbosity = val; }
		inline void set_num_threads(int num) { num_threads = num; }
		// Simple getters
		inline bool get_recover() const { return recover; }
		inline unsigned long int get_random_seed() const { return rngseed; }
		inline std::string const& get_datafile_name() const { return datafile_name; }
		inline std::string const& get_wisdom_file_name() const { return wisdom_file_name; }
		inline std::string const& get_copy_from() const { return copy_from; }
		inline SaveWhat get_save_what() const { return save_what; }
		inline bool get_clobber() const { return clobber; }
		inline int get_verbosity() const { return verbosity; }
		inline size_t get_num_threads() const { return num_threads; }
		inline size_t get_sizex() const { return sizex; }
		inline size_t get_sizey() const { return sizey; }
		inline double get_lenx() const { return lenx; }
		inline BoundaryType get_boundary_type() const { return boundary; }
		inline double get_grid_delta() const { return lenx/static_cast<double>(sizex); }
		inline size_t get_N() const { return N; }
		inline OrthoAlgorithm get_ortho_algorithm() const { return ortho_alg; }
		inline unsigned int get_fftw_flags() const { return fftw_flags; }
		inline InitialStatePreset get_initialstate_preset () const { return initialstate_preset; }
		inline initialstatefunc get_initialstate_func() const { return initialstate_func; }
		inline std::string const& get_initialstate_description() const { return initialstate_description; }
		inline std::string const& get_potential_type() const { return potential_type; }
		inline ConvergenceTest const& get_timestep_convergence_test() const { return *timestep_convergence_test; }
		inline ConvergenceTest const& get_final_convergence_test() const { return *final_convergence_test; }
		inline std::string const& get_noise_type() const { return noise_type; }
		inline std::string const& get_noise_constraint_type() const { return noise_constraint_type; }
		inline double get_B() const { return B; }
		inline int get_halforder() const { return halforder; }
		inline std::list<double> const& get_eps_values() const { return eps_values; }
		inline double get_eps_divisor() const { return eps_divisor; }
		inline bool get_exhaust_eps() const { return exhaust_eps; }
		inline size_t get_needed_to_converge() const { return needed_to_converge; }
		inline size_t get_ignore_lowest() const { return ignore_lowest; }
		inline int get_max_steps() const { return max_steps; }
		inline double get_min_time_step() const { return min_time_step; }
		/*
		 * Default values
		 */
		static const bool default_recover;
		static const unsigned long int default_rngseed;
		static const char default_datafile_name[];
		static const char default_wisdom_file_name[];
		static const SaveWhat default_save_what;
		static const bool default_clobber;
		static const int default_verbosity;
		static const size_t default_num_threads;
		static const BoundaryType default_boundary;
		static const size_t default_sizex;
		static const size_t default_sizey;
		static const double default_lenx;
		static const size_t default_N;
		static const OrthoAlgorithm default_ortho_alg;
		static const InitialStatePreset default_initialstate_preset;
		static const char default_potential_type[];
		static const char default_timestep_convergence_test_string[];
		static const char default_final_convergence_test_string[];
		static const char default_noise_type[];
		static const char default_noise_constraint_type[];
		static const double default_B;
		static const int default_halforder;
		static const double default_initial_eps;
		static const double default_eps_divisor;
		static const bool default_exhaust_eps;
		static const size_t default_needed_to_converge;
		static const size_t default_ignore_lowest;
		static const int default_max_steps;
		static const double default_min_time_step;
	private:
		void set_to_defaults();
		/* 
		 * Data members
		 */
		bool recover;
		// RNG seed
		unsigned long int rngseed;
		// Datafile parameters
		std::string datafile_name;
		std::string wisdom_file_name;
		std::string copy_from;
		SaveWhat save_what;
		bool clobber;	// If true, an existing datafile will be overwritten
		// Output parameters
		int verbosity;
		// General performance parameters
		size_t num_threads;
		OrthoAlgorithm ortho_alg;
		unsigned int fftw_flags;
		// Grid parameters
		BoundaryType boundary;
		size_t sizex;
		size_t sizey;
		double lenx;
		// State parameters
		size_t N;
		InitialStatePreset initialstate_preset;
		std::string initialstate_description;
		initialstatefunc initialstate_func;
		// Potential & field parameters
		std::string potential_type;
		double B;
		std::string noise_type;
		std::string noise_constraint_type;
		// Time evolution parameters
		int halforder;					// Half the order of operator factorization
		std::list<double> eps_values;	// Time step values to use
		double eps_divisor;				// Divisor used to generate more time step values once eps_values is exhausted
		bool exhaust_eps;				// If true, run just one iteration with each value of eps_values and then continue normally
		// Convergence criteria
		ConvergenceTest const* timestep_convergence_test;
		ConvergenceTest const* final_convergence_test;
		size_t needed_to_converge;
		size_t ignore_lowest;
		// Bailout criteria
		int max_steps;
		double min_time_step;
};

inline void Parameters::set_timestep_convergence_test(std::string const& str) {
	delete timestep_convergence_test;
	timestep_convergence_test = parse_convergence_description(str);
}

inline void Parameters::set_final_convergence_test(std::string const& str) {
	delete final_convergence_test;
	final_convergence_test = parse_convergence_description(str);
}

std::ostream& operator<<(std::ostream& stream, const Parameters& params);

#endif // _PARAMETERS_HPP_
