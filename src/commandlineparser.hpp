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
 * A command line parsing class for itp2d.
 * This is essentially a wrapper over the several argument classes provided by
 * TCLAP, the Templatized C++ Command Line Parser Library
 */

#ifndef _COMMANDLINEPARSER_HPP_
#define _COMMANDLINEPARSER_HPP_

#include <tclap/CmdLine.h>

#include "itp2d_common.hpp"
#include "exceptions.hpp"
#include "parameters.hpp"
#include "rng.hpp"

class CommandLineParser {
	public:
		CommandLineParser();										// The basic usage:
		void parse(std::vector<std::string>& args);					// parse a string of command line arguments
		Parameters const& get_params() const { return params; }		// ...and return the corresponding Parameters instance
		// Documentation strings for each command line parameter. These are printed with --help
		static const char help_highmem[];
		static const char help_wisdom_file_name[];
		static const char help_noise[];
		static const char help_noise_constraint[];
		static const char help_recover[];
		static const char help_rngseed[];
		static const char help_min_time_step[];
		static const char help_max_steps[];
		static const char help_exhaust_eps_values[];
		static const char help_eps_divisor[];
		static const char help_eps_values[];
		static const char help_order[];
		static const char help_N[];
		static const char help_dirichlet[];
		static const char help_pi[];
		static const char help_lenx[];
		static const char help_sizey[];
		static const char help_sizex[];
		static const char help_size[];
		static const char help_B[];
		static const char help_ignore_lowest[];
		static const char help_needed_to_converge[];
		static const char help_num_threads[];
		static const char help_quietness[];
		static const char help_verbosity[];
		static const char help_save_nothing[];
		static const char help_save_onlyenergies[];
		static const char help_save_everything[];
		static const char help_clobber[];
		static const char help_copy_from[];
		static const char help_datafile_name[];
		static const char help_timestep_convtest[];
		static const char help_final_convtest[];
		static const char help_potential[];
		// Epilogue to be printed after the command line parameter descriptions
		static const char help_epilogue[];
	private:
		template <typename Type> inline void throw_if_nonpositive(Type value, std::string const& name);
		template <typename Type> inline void throw_if_nonpositive(TCLAP::ValueArg<Type>& arg);
		template <typename Type> inline void throw_if_negative(TCLAP::ValueArg<Type>& arg);
		Parameters params;
		TCLAP::CmdLine cmd;
		TCLAP::SwitchArg arg_highmem;
		TCLAP::ValueArg<std::string> arg_wisdom_file_name;
		TCLAP::ValueArg<std::string> arg_noise;
		TCLAP::ValueArg<std::string> arg_noise_constraint;
		TCLAP::SwitchArg arg_recover;
		TCLAP::ValueArg<unsigned long int> arg_rngseed;
		TCLAP::ValueArg<double> arg_min_time_step;
		TCLAP::ValueArg<int> arg_max_steps;
		TCLAP::SwitchArg arg_exhaust_eps_values;
		TCLAP::ValueArg<double> arg_eps_divisor;
		TCLAP::MultiArg<double> arg_eps_values;
		TCLAP::ValueArg<int> arg_order;
		TCLAP::ValueArg<size_t> arg_N;
		TCLAP::SwitchArg arg_dirichlet;
		TCLAP::SwitchArg arg_pi;
		TCLAP::ValueArg<double> arg_lenx;
		TCLAP::ValueArg<size_t> arg_sizey;
		TCLAP::ValueArg<size_t> arg_sizex;
		TCLAP::ValueArg<size_t> arg_size;
		TCLAP::ValueArg<double> arg_B;
		TCLAP::ValueArg<size_t> arg_ignore_lowest;
		TCLAP::ValueArg<size_t> arg_needed_to_converge;
		TCLAP::ValueArg<size_t> arg_num_threads;
		TCLAP::MultiSwitchArg arg_quietness;
		TCLAP::MultiSwitchArg arg_verbosity;
		TCLAP::SwitchArg arg_save_nothing;
		TCLAP::SwitchArg arg_save_onlyenergies;
		TCLAP::SwitchArg arg_save_everything;
		TCLAP::SwitchArg arg_clobber;
		TCLAP::ValueArg<std::string> arg_copy_from;
		TCLAP::ValueArg<std::string> arg_datafile_name;
		TCLAP::ValueArg<std::string> arg_timestep_convtest;
		TCLAP::ValueArg<std::string> arg_final_convtest;
		TCLAP::ValueArg<std::string> arg_potential;
};

/* Bunch of helper functions for throwing exceptions for invalid command line
 * arguments
 */

template <typename Type>
inline void CommandLineParser::throw_if_nonpositive(Type value, std::string const& name) {
	if (value <= 0)
		throw TCLAP::CmdLineParseException("Non-positive value not allowed.", name);
}

template <typename Type>
inline void CommandLineParser::throw_if_nonpositive(TCLAP::ValueArg<Type>& arg) {
	throw_if_nonpositive(arg.getValue(), arg.getName());
}

template <typename Type>
inline void CommandLineParser::throw_if_negative(TCLAP::ValueArg<Type>& arg) {
	if (arg.getValue() < 0)
		throw TCLAP::CmdLineParseException("Negative value not allowed.", arg.getName());
}

#endif // _COMMANDLINEPARSER_HPP_
