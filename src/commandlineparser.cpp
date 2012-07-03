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

// Documentation strings

const char CommandLineParser::help_highmem[] = "\
Use a different orthonormalization algorithm, which doubles the memory usage but *possibly* offers \
better performance.";

const char CommandLineParser::help_wisdom_file_name[] = "\
File name to use for FFTW wisdom.";

const char CommandLineParser::help_noise[] = "\
Description of possible noise added to the potential. Valid descriptions:\n\
Gaussian spikes with the prescribed density and normally distributed amplitude and width:\n\
\tgaussian(density,amp_mean,width_mean)\n\
\tgaussian(density,amp_mean,amp_stdev,width_mean,width_stdev)\n\
See header noise.hpp for details.";

const char CommandLineParser::help_recover[] = "\
Restart simulation instead of quitting on some fatal errors.";

const char CommandLineParser::help_rngseed[] = "\
Provide a seed for the random number generator. If not set, one is generated based on the current time.";

const char CommandLineParser::help_min_time_step[] = "\
Bail out if the imaginary time step goes below this value.";

const char CommandLineParser::help_max_steps[] = "\
Bail out after this many iterations.";

const char CommandLineParser::help_exhaust_eps_values[] = "\
Exhaust user-specified (with argument --timestep) list of time step values before starting \
convergence checking, i.e. do one iteration with each provided time step value.";

const char CommandLineParser::help_eps_divisor[] = "\
Divisor used to decrease the time step once all user-specified values have been used.";

const char CommandLineParser::help_eps_values[] = "\
A value for the initial imaginary time step. You can give this argument multiple times to specify a \
list of values that will be used in the order you specify them.";

const char CommandLineParser::help_order[] = "\
Order of operator splitting. Has to be an even number.";

const char CommandLineParser::help_N[] = "\
Number of states to use in computations. If not set, defaults to the value of --states plus \
additional 25%.";

const char CommandLineParser::help_dirichlet[] = "\
Use Dirichlet boundary conditions. By default itp2d uses periodic boundary conditions.";

const char CommandLineParser::help_pi[] = "\
Multiply grid length given by --lenx by pi, i.e., use '--pi -l 1.0' to get a pi by pi box.";

const char CommandLineParser::help_lenx[] = "\
Grid length in x-direction.";

const char CommandLineParser::help_sizey[] = "\
Number of grid points in the y-direction.";

const char CommandLineParser::help_sizex[] = "\
Number of grid points in the x-direction.";

const char CommandLineParser::help_size[] = "\
Number of grid points along each dimension. Sets the same value for x- and y-directions.";

const char CommandLineParser::help_B[] = "\
Strength of the external magnetic field.";

const char CommandLineParser::help_ignore_lowest[] = "\
Ignore this many lowest states in convergence checking.";

const char CommandLineParser::help_needed_to_converge[] = "\
Number of states wanted to converge.";

const char CommandLineParser::help_num_threads[] = "\
Use this many threads.";

const char CommandLineParser::help_quietness[] = "\
Decrease verbosity of output.";

const char CommandLineParser::help_verbosity[] = "\
Increase verbosity of output.";

const char CommandLineParser::help_save_nothing[] = "\
Run without saving anything on disk.";

const char CommandLineParser::help_save_onlyenergies[] = "\
Save only final state energies, not the states themselves.";

const char CommandLineParser::help_save_everything[] = "\
Save state data after each step. Causes MASSIVE datafiles.";

const char CommandLineParser::help_clobber[] = "\
Overwrite datafile if it exists.";

const char CommandLineParser::help_copy_from[] = "\
Copy state data from specified datafile.";

const char CommandLineParser::help_datafile_name[] = "\
File name to save data to.";

const char CommandLineParser::help_timestep_convtest[] = "\
Description of the test for timestep convergence. Valid descriptions:\n\
No convergence checking:\n\
\tnone\n\
One-step timestep convergence:\n\
\tonestep\n\
Absolute energy change less than value:\n\
\tabsEchange(value)\n\
Relative energy change less than value:\n\
\trelEchange(value)\n\
Standard deviation of energy less than value or change of standard \
deviation less than value2:\n\
\tabsdeviation(value)\n\
\tabsdeviation(value,value2)\n\
Relative standard deviation of energy less than value or change of relative \
standard deviation less than value2:\n\
\tdeviation(value)\n\
\tdeviation(value,value2)\n\
In all cases 'change' means change between successive iterations.\n\
See header convergence.hpp for details.";

const char CommandLineParser::help_final_convtest[] = "\
Description of the test for final convergence. \
See the documentation for --timestep_convtest for details.";

const char CommandLineParser::help_potential[] = "\
Description of the potential. Valid descriptions:\n\
Zero potential:\n\
\tzero\n\
Harmonic oscillator with frequency w, centered at (x0,y0):\n\
\tharmonic(w)\n\
\tharmonic(w,x0,y0)\n\
Square box with power function walls:\n\
\tprettyhardsquare(exponent)\n\
Soft-walled pentagon:\n\
\tsoftpentagon\n\
The Henon-Heiles potential:\n\
\thenonheiles(a,b)\n\
Gaussian blob:\n\
\tgaussian(amplitude,width)\n\
\tgaussian(amplitude,width,x0,y0)\n\
Quartic oscillator potential (x^2 * y^2)/2 + b(x^4 + y^4)/4, rotated by pi/4:\n\
\tquartic(b)\n\
See header potential.hpp for details.";

const char CommandLineParser::help_epilogue[] = "\
All values given or received by itp2d are in SI-based Hartree atomic units.\n\
Copyright 2012 Perttu Luukko\n\
itp2d is free software: you can redistribute it and/or modify it under the terms of the GNU General \
Public License as published by the Free Software Foundation, either version 3 of the License, or (at \
your option) any later version. itp2d is distributed in the hope that it will be useful, but \
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A \
PARTICULAR PURPOSE. See the GNU General Public License for more details.";

// The constructor. Simply initialize all the argument instances with the help strings and default values.

CommandLineParser::CommandLineParser() :
	params(),
	cmd(help_epilogue, ' ', version_string),
	arg_highmem("", "highmem-orthonormalization", help_highmem, cmd),
	arg_wisdom_file_name("", "wisdomfile", help_wisdom_file_name, false, Parameters::default_wisdom_file_name, "FILENAME", cmd),
	arg_noise("", "noise", help_noise, false, Parameters::default_noise_type, "STRING", cmd),
	arg_recover("", "recover", help_recover, cmd),
	arg_rngseed("", "rngseed", help_rngseed, false, Parameters::default_rngseed, "NUM", cmd),
	arg_min_time_step("", "mineps", help_min_time_step, false, Parameters::default_min_time_step, "FLOAT", cmd),
	arg_max_steps("", "maxsteps", help_max_steps, false, Parameters::default_max_steps, "NUM", cmd),
	arg_exhaust_eps_values("", "exhaust-timestep-list", help_exhaust_eps_values, cmd),
	arg_eps_divisor("D", "timestep-divisor", help_eps_divisor, false, Parameters::default_eps_divisor, "FLOAT", cmd),
	arg_eps_values("e", "timestep", help_eps_values, false, "FLOAT", cmd),
	arg_order("d", "order", help_order, false, 2*Parameters::default_halforder, "NUM", cmd),
	arg_N("N", "totalstates", help_N, false, Parameters::default_N, "NUM", cmd),
	arg_dirichlet("", "dirichlet", help_dirichlet, cmd),
	arg_pi("", "pi", help_pi, cmd),
	arg_lenx("l", "lenx", help_lenx, false, Parameters::default_lenx, "FLOAT", cmd),
	arg_sizey("y", "sizey", help_sizey, false, Parameters::default_sizey, "NUM", cmd),
	arg_sizex("x", "sizex", help_sizex, false, Parameters::default_sizex, "NUM", cmd),
	arg_size("s", "size", help_size, false, Parameters::default_sizex, "NUM", cmd),
	arg_B("B", "magneticfield", help_B, false, Parameters::default_B, "FLOAT", cmd),
	arg_ignore_lowest("", "ignore-lowest", help_ignore_lowest, false, Parameters::default_ignore_lowest, "NUM", cmd),
	arg_needed_to_converge("n", "states", help_needed_to_converge, false, Parameters::default_needed_to_converge, "NUM", cmd),
	arg_num_threads("t", "threads", help_num_threads, false, Parameters::default_num_threads, "NUM", cmd),
	arg_quietness("q", "quiet", help_quietness, cmd),
	arg_verbosity("v", "verbose", help_verbosity, cmd),
	arg_save_nothing("", "save-nothing", help_save_nothing, cmd),
	arg_save_onlyenergies("", "save-only-energies", help_save_onlyenergies, cmd),
	arg_save_everything("", "save-everything", help_save_everything, cmd),
	arg_clobber("f", "force", help_clobber, cmd),
	arg_copy_from("", "copy-states", help_copy_from, false, "", "FILENAME", cmd),
	arg_datafile_name("o", "datafile", help_datafile_name, false, Parameters::default_datafile_name, "FILENAME", cmd),
	arg_timestep_convtest("T", "timestep-convergence-test", help_timestep_convtest, false, Parameters::default_timestep_convergence_test_string, "STRING", cmd),
	arg_final_convtest("F", "final-convergence-test", help_final_convtest, false, Parameters::default_final_convergence_test_string, "STRING", cmd),
	arg_potential("p", "potential", help_potential, false, Parameters::default_potential_type, "STRING", cmd) {}

void CommandLineParser::parse(std::vector<std::string>& args) {
	cmd.setExceptionHandling(false);
	// Run command line parser
	cmd.parse(args);
	// Do some validations on the arguments given
	{
		int save_args_counter = 0;
		if (arg_save_nothing.isSet())
			save_args_counter++;
		if (arg_save_onlyenergies.isSet())
			save_args_counter++;
		if (arg_save_everything.isSet())
			save_args_counter++;
	if (save_args_counter > 1)
		throw TCLAP::CmdLineParseException("Arguments cannot be set together.",
				arg_save_nothing.getName()+", "+
				arg_save_onlyenergies.getName()+", "+
				arg_save_everything.getName());
	}
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
	if (arg_rngseed.isSet())
		params.set_random_seed(arg_rngseed.getValue());
	else
		params.set_random_seed(RNG::produce_random_seed());
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
	if (arg_save_nothing.getValue())
		params.save_what = Parameters::Nothing;
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
