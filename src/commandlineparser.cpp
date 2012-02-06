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

#include "commandlineparser.hpp"

CommandLineParser::CommandLineParser() :
	params(),
	cmd("All values given or received by itp2d are in SI-based Hartree atomic units.\nCopyright 2012 Perttu Luukko\nitp2d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. itp2d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.", ' ', version_string),
	arg_highmem("", "highmem", "Use a different orthonormalization algorithm, which doubles the memory usage but *possibly* offers better performance.", cmd),
	arg_wisdom_file_name("", "wisdom", "File name to use for FFTW wisdom.", false, Parameters::default_wisdom_file_name, "FILENAME", cmd),
	arg_noise("", "noise", "Description of possible noise added to the potential. Valid descriptions:\nGaussian spikes with the prescribed density and normally distributed amplitude and width:\n\tgaussians(density,amp_mean,amp_stdev,width_mean,width_stdev)\nSee header noise.hpp for details.", false, Parameters::default_noise_type, "STRING", cmd),
	arg_recover("", "recover", "Restart simulation instead of quitting on some fatal errors.", cmd),
	arg_min_time_step("", "mineps", "Bail out if the imaginary time step goes below this value.", false, Parameters::default_min_time_step, "FLOAT", cmd),
	arg_max_steps("", "maxsteps", "Bail out after this many iterations.", false, Parameters::default_max_steps, "NUM", cmd),
	arg_exhaust_eps_values("", "exhaust_eps", "Exhaust user-specified (with argument --timestep) list of time step values before starting convergence checking, i.e. do one iteration with each provided time step value.", cmd),
	arg_eps_divisor("D", "eps_divisor", "Divisor used to decrease the time step once all user-specified values have been used.", false, Parameters::default_eps_divisor, "FLOAT", cmd),
	arg_eps_values("e", "timestep", "A value for the initial imaginary time step. You can give this argument multiple times to specify a list of values that will be used in the order you specify them.", false, "FLOAT", cmd),
	arg_order("d", "order", "Order of operator splitting. Has to be an even number.", false, 2*Parameters::default_halforder, "NUM", cmd),
	arg_N("N", "totalstates", "Number of states to use in computations. If not set, defaults to the value of --states plus additional 25%.", false, Parameters::default_N, "NUM", cmd),
	arg_dirichlet("", "dirichlet", "Use Dirichlet boundary conditions. By default itp2d uses periodic boundary conditions.", cmd),
	arg_pi("", "pi", "Multiply grid length given by --lenx by pi, i.e., use '--pi -l 1.0' to get a pi by pi box.", cmd),
	arg_lenx("l", "lenx", "Grid length in x-direction.", false, Parameters::default_lenx, "FLOAT", cmd),
	arg_sizey("y", "sizey", "Number of grid points in the y-direction.", false, Parameters::default_sizey, "NUM", cmd),
	arg_sizex("x", "sizex", "Number of grid points in the x-direction.", false, Parameters::default_sizex, "NUM", cmd),
	arg_size("s", "size", "Number of grid points along each dimension. Sets the same value for x- and y-directions.", false, Parameters::default_sizex, "NUM", cmd),
	arg_B("B", "magneticfield", "Strength of the external magnetic field.", false, Parameters::default_B, "FLOAT", cmd),
	arg_ignore_lowest("", "ignorelow", "Ignore this many lowest states in convergence checking.", false, Parameters::default_ignore_lowest, "NUM", cmd),
	arg_needed_to_converge("n", "states", "Number of states wanted to converge.", false, Parameters::default_needed_to_converge, "NUM", cmd),
	arg_num_threads("t", "threads", "Use this many threads.", false, Parameters::default_num_threads, "NUM", cmd),
	arg_quietness("q", "quiet", "Decrease verbosity of output.", cmd),
	arg_verbosity("v", "verbose", "Increase verbosity of output.", cmd),
	arg_save_onlyenergies("", "onlyenergies", "Save only final state energies, not the states themselves.", cmd),
	arg_save_everything("", "everything", "Save state data after each step. Causes MASSIVE datafiles.", cmd),
	arg_clobber("f", "force", "Overwrite datafile if it exists.", cmd),
	arg_copy_from("", "copystates", "Copy state data from specified datafile.", false, "", "FILENAME", cmd),
	arg_datafile_name("o", "datafile", "File name to save data to.", false, Parameters::default_datafile_name, "FILENAME", cmd),
	arg_timestep_convtest("T", "timestep_convtest", "Description of the test for timestep convergence. Valid descriptions:\nNo convergence checking:\n\tnone\nOne-step timestep convergence:\n\tonestep\nAbsolute energy change less than value:\n\tabsEchange(value)\nRelative energy change less than value:\n\trelEchange(value)\nRelative standard deviation of energy less than value or change of relative standard deviation less than value2:\n\tdeviation(value,value2)\nIn all cases 'change' means change between successive iterations.\nSee header convergence.hpp for details.", false, Parameters::default_timestep_convergence_test_string, "STRING", cmd),
	arg_final_convtest("F", "final_convtest", "Description of the test for final convergence. See the documentation for --timestep_convtest for details.", false, Parameters::default_final_convergence_test_string, "STRING", cmd),
	arg_potential("p", "potential", "Description of the potential. Valid descriptions:\nZero potential:\n\tzero\nHarmonic oscillator with frequency w:\n\tharmonic(w)\nSquare box with power function walls:\n\tprettyhardsquare(exponent)\nSoft-walled pentagon:\n\tsoftpentagon\nThe Henon-Heiles potential:\n\thenonheiles(a,b)\nGaussian blob:\n\tgaussian(amplitude,width,x0,y0)\nSee header potential.hpp for details.", false, "default", "STRING", cmd) {}

void CommandLineParser::parse(std::vector<std::string>& args) {
	try {
		// Run command line parser
		cmd.parse(args);
		// Do some validations on the arguments given
		if (arg_save_everything.isSet() and arg_save_onlyenergies.isSet())
			throw TCLAP::CmdLineParseException("Both arguments cannot be set together.", arg_save_everything.getName()+" and "+arg_save_onlyenergies.getName());
		throw_if_nonpositive(arg_num_threads);
		if (arg_size.isSet() and (arg_sizex.isSet() or arg_sizey.isSet()))
			throw TCLAP::CmdLineParseException("Arguments cannot be set together.",
				arg_size.getName()+" and ("+arg_sizey.getName()+" or "+arg_sizex.getName()+")");
		throw_if_nonpositive(arg_size);
		throw_if_nonpositive(arg_sizex);
		throw_if_nonpositive(arg_sizey);
		throw_if_nonpositive(arg_lenx);
		std::vector<double> const& eps_values = arg_eps_values.getValue();
		for (std::vector<double>::const_iterator it = eps_values.begin(); it != eps_values.end(); it++) {
			throw_if_nonpositive(*it, arg_eps_values.getName());
		}
		throw_if_nonpositive(arg_eps_divisor);
		throw_if_nonpositive(arg_N);
		throw_if_nonpositive(arg_needed_to_converge);
		if (arg_N.isSet() and arg_needed_to_converge.isSet() and arg_N.getValue() < (arg_ignore_lowest.getValue() + arg_needed_to_converge.getValue()))
			throw TCLAP::CmdLineParseException("Number of states has to be large enough to include at least all the states wanted to converge!", arg_N.getName());
		throw_if_nonpositive(arg_order);
		if (arg_order.getValue() % 2 != 0)
			throw TCLAP::CmdLineParseException("Has to be even.", arg_order.getName());
		throw_if_negative(arg_min_time_step);
		throw_if_nonpositive(arg_max_steps);
		// Build the Parameters class instance based on the command line options given
		params.recover = arg_recover.getValue();
		params.datafile_name = arg_datafile_name.getValue();
		params.set_wisdom_file_name(arg_wisdom_file_name.getValue());
		params.copy_from = arg_copy_from.getValue();
		if (arg_copy_from.isSet())
			params.initialstate_preset = Parameters::CopyFromFile;
		params.potential_type = arg_potential.getValue();
		params.set_timestep_convergence_test(arg_timestep_convtest.getValue());
		params.set_final_convergence_test(arg_final_convtest.getValue());
		params.set_noise(arg_noise.getValue());
		if (arg_save_everything.getValue())
			params.save_what = Parameters::Everything;
		if (arg_save_onlyenergies.getValue())
			params.save_what = Parameters::OnlyEnergies;
		params.clobber = arg_clobber.getValue();
		params.verbosity = Parameters::default_verbosity + arg_verbosity.getValue() - arg_quietness.getValue();
		params.num_threads = arg_num_threads.getValue();
		params.sizex = arg_sizex.getValue();
		params.sizey = arg_sizey.getValue();
		if (arg_size.isSet()) {
			params.sizex = arg_size.getValue();
			params.sizey = arg_size.getValue();
		}
		params.lenx = ((arg_pi.isSet())? pi : 1.0) *arg_lenx.getValue();
		params.boundary = arg_dirichlet.getValue()? Dirichlet : Periodic;
		params.ortho_alg = arg_highmem.getValue()? HighMem : Default;
		for (std::vector<double>::const_iterator it = eps_values.begin(); it != eps_values.end(); it++) {
			params.add_eps_value(*it);
		}
		params.eps_divisor = arg_eps_divisor.getValue();
		params.exhaust_eps = arg_exhaust_eps_values.getValue();
		if (not arg_needed_to_converge.isSet() and arg_N.isSet())
			params.needed_to_converge = arg_N.getValue();
		else
			params.needed_to_converge = arg_needed_to_converge.getValue();
		// If the total number of states is not given, set it to the number required to converge plus 25%
		if (arg_needed_to_converge.isSet() and not arg_N.isSet())
			params.N = static_cast<size_t>(static_cast<double>(arg_needed_to_converge.getValue())*1.25);
		else
			params.N = arg_N.getValue();
		params.ignore_lowest = arg_ignore_lowest.getValue();
		params.B = arg_B.getValue();
		params.halforder = arg_order.getValue()/2;
		params.min_time_step = arg_min_time_step.getValue();
		params.max_steps = arg_max_steps.getValue();
	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "Command line parsing error: " << e.error() << " " << e.argId() << std::endl;
		throw GeneralError("Command line parsing failed.");
	}
}
