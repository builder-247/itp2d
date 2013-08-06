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
 * A complete simulation program using the imaginary time propagation algorithm
 * to solve the time-independent Schrödinger equation in 2D. This basically
 * just wraps together the ITPSystem class with some command line parameter
 * parsing from CommandLineParser.
 */

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <csignal>

#include <unistd.h>

#include "parameters.hpp"
#include "commandlineparser.hpp"
#include "itpsystem.hpp"

using namespace std;

// Global flags for signal handling
volatile sig_atomic_t abort_flag = false;
volatile sig_atomic_t save_flag = false;

const char abort_flag_note[] = "\nCaught SIGINT. Saving data and quitting at next convenient spot.\nPress Ctrl-C again to signal immediate stop.\n";
const size_t abort_flag_note_size = sizeof(abort_flag_note)+1;
const char save_flag_note[] = "\nCaught SIGUSR1. Saving states at next convenient stop and continuing.\n";
const size_t save_flag_note_size = sizeof(save_flag_note)+1;

// Signal handlers for SIGINT and SIGUSR1
void sigint_handler(__attribute__((unused)) int s) {
	if (abort_flag)
		exit(1);
	write(STDERR_FILENO, abort_flag_note, abort_flag_note_size);
	abort_flag = true;
}

void sigusr1_handler(__attribute__((unused)) int s) {
	write(STDERR_FILENO, save_flag_note, save_flag_note_size);
	save_flag = true;
}

int main(int argc, char* argv[]) {
	// Trap SIGINT and SIGUSR1
	signal(SIGINT, sigint_handler);
	signal(SIGUSR1, sigusr1_handler);
	// Parse parameters
	vector<string> args(argv, argv+argc);
	string program_name(args.front()); // TCLAP eats the first element of args
	CommandLineParser parser;
	try {
		parser.parse(args);
	}
	catch (TCLAP::ArgException &e) {
		cerr << "Command line parsing error:" << endl
			<< "\tError: " << e.error() << endl
			<< "\t" << e.argId() << endl << endl
			<< "For documentation on what command line arguments are available" << endl
			<< "and what they mean, please type:" << endl
			<< program_name << " --help" << endl;
		return 2;
	}
	catch (TCLAP::ExitException &e) {
		return e.getExitStatus();
	}
	Parameters params(parser.get_params());
	// Import FFTW Wisdom if available
	std::string const& fftw_wisdom_filename = params.get_wisdom_file_name();
	FILE* wisdom_file = fopen(fftw_wisdom_filename.c_str(), "r");
	if (wisdom_file != NULL) {
		fftw_import_wisdom_from_file(wisdom_file);
		fclose(wisdom_file);
	}
	// Initialize ITPSystem
	ITPSystem* sys = new ITPSystem(params, &abort_flag, &save_flag);
	// Main loop
	while(not sys->is_finished()) {
		sys->step();
	}
	// Save FFTW Wisdom
	wisdom_file = fopen(fftw_wisdom_filename.c_str(), "w");
	fftw_export_wisdom_to_file(wisdom_file);
	// Cleanup and exit
	const bool error_flag = sys->get_error_flag();
	delete sys;
	fftw_cleanup();
	if (error_flag)
		return 1;
	else
		return 0;
}
