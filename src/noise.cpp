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

#include "noise.hpp"

// default parameters

const double GaussianNoise::default_relative_amplitude_stdev = 0.20;
const double GaussianNoise::default_relative_width_stdev = 0.20;

const double CoulombImpurities::default_exponent = 1.0;
const double CoulombImpurities::default_alpha = 1.0;
const double CoulombImpurities::default_maxd = 1.0;

// noise parser function

Noise const* parse_noise_description(std::string const& str) {
	// Parse with the generic parse_parameter_string function and extract the
	// name and parameters.
	name_parameters_pair p;
	try {
		p = parse_parameter_string(str);
	}
	catch (ParseError& e) {
		std::cerr << e.what() << std::endl;
		throw InvalidNoiseType(str);
	}
	std::string const& name = p.first;
	std::vector<double> const& params = p.second;
	// Simply delegate to the individual constructors based on name
	if (name == "no" or name == "none" or name == "zero")
		return new NoNoise(params);
	else if (name == "gaussian" or name == "gaussians" or name == "gaussiannoise")
		return new GaussianNoise(params);
	else if (name == "coulomb" or name == "coulombimpurities")
		return new CoulombImpurities(params);
	else
		throw UnknownNoiseType(str);
}

// NoNoise

NoNoise::NoNoise(std::vector<double> params) {
	if (not params.empty())
		throw InvalidNoiseType("Noise type NoNoise does not take parameters");
	init();
}

// GaussianNoise

GaussianNoise::GaussianNoise(std::vector<double> params) {
	if (params.size() == 3) {
		density = params[0];
		amplitude_mean = params[1];
		width_mean = params[2];
		amplitude_stdev = default_relative_amplitude_stdev * amplitude_mean;
		width_stdev = default_relative_width_stdev * width_mean;
	}
	else if (params.size() == 5) {
		density = params[0];
		amplitude_mean = params[1];
		amplitude_stdev = params[2];
		width_mean = params[3];
		width_stdev = params[4];
	}
	else
		throw InvalidNoiseType("Noise type GaussianNoise takes either 3 or 5 parameters");
	init();
}

void GaussianNoise::add_noise(DataLayout const& dl, double* pot_values, RNG& rng, Constraint const& constraint) const {
	// The number of impurities on the calculation plane is Poisson distributed around the mean value
	const double lambda = density*dl.lenx*dl.leny;
	const unsigned int N = rng.poisson_rand(lambda);
	// First, randomly distribute the centers, amplitudes and widths of the spikes
	typedef std::tr1::tuple<double, double, double, double> spike; // x-coord, y-coord, amplitude and width of a spike
	typedef std::vector<spike> spike_vector;
	spike_vector spikes;
	for (unsigned int n=0; n<N; n++) {
		// First randomize the position
		const double x = (rng.uniform_rand()-0.5)*dl.lenx;
		const double y = (rng.uniform_rand()-0.5)*dl.leny;
		if (not constraint.check(x, y)) {
			continue;
		}
		// Determine amplitude and width
		const double A = amplitude_mean + amplitude_stdev*rng.gaussian_rand();
		const double w = width_mean + width_stdev*rng.gaussian_rand();
		spikes.push_back(spike(x, y, A, w));
	}
	// Then, loop through all spikes and add its effect into pot_values. This
	// is not a very efficient way to do this, but this only needs to be done
	// once so it doesn't matter.
	for (spike_vector::const_iterator it = spikes.begin(); it != spikes.end(); ++it) {
		const double sx = std::tr1::get<0>(*it);
		const double sy = std::tr1::get<1>(*it);
		const double A = std::tr1::get<2>(*it);
		const double w = std::tr1::get<3>(*it);
		const double w2 = w*w;
		for (size_t x=0; x<dl.sizex; x++) {
			const double px = dl.get_posx(x);
			for (size_t y=0; y<dl.sizey; y++) {
				const double py = dl.get_posy(y);
				const double rx = sx-px;
				const double ry = sy-py;
				const double r2 = rx*rx + ry*ry;
				dl.value(pot_values, x, y) += A*exp(-0.5*(r2/w2));
			}
		}
	}
}

void GaussianNoise::init() {
	std::stringstream ss;
	ss << "gaussian spikes, density = " << density << ", amplitude ~ N(" << amplitude_mean
		<< "," << amplitude_stdev << "^2), width ~ N(" << width_mean << "," << width_stdev << "^2)";
	description = ss.str();
}

// CoulombImpurities

CoulombImpurities::CoulombImpurities(std::vector<double> params) {
	// It looks like this "design pattern" has reached its end
	if (params.size() == 4) {
		density = params[0];
		exponent = params[1];
		alpha = params[2];
		maxd = params[3];
	}
	else if (params.size() == 3) {
		density = params[0];
		exponent = params[1];
		alpha = params[2];
		maxd = default_maxd;
	}
	else if (params.size() == 2) {
		density = params[0];
		exponent = params[1];
		alpha = default_alpha;
		maxd = default_maxd;
	}
	else if (params.size() == 1) {
		density = params[0];
		exponent = default_exponent;
		alpha = default_alpha;
		maxd = default_maxd;
	}
	else
		throw InvalidNoiseType("Noise type CoulombImpurities takes 1-4 parameters");
	init();
}

void CoulombImpurities::add_noise(DataLayout const& dl, double* pot_values, RNG& rng, Constraint const& constraint) const {
	// The number of impurities inside the box is Poisson distributed around the mean value
	const double lambda = density*dl.lenx*dl.leny*2*maxd;
	const unsigned int N = rng.poisson_rand(lambda);
	// First, randomly distribute the positions of the impurities
	typedef std::tr1::tuple<double, double, double> coords; // x-coord, y-coord, z-coord
	typedef std::vector<coords> coords_vector;
	coords_vector impurity_coords;
	for (unsigned int n=0; n<N; n++) {
		const double x = (rng.uniform_rand()-0.5)*dl.lenx;
		const double y = (rng.uniform_rand()-0.5)*dl.leny;
		if (not constraint.check(x, y)) {
			continue;
		}
		const double z = (rng.uniform_rand()-0.5)*2*maxd;
		impurity_coords.push_back(coords(x, y, z));
	}
	// Then, loop through all impurities and add its effect into pot_values. This
	// is not a very efficient way to do this, but this only needs to be done
	// once so it doesn't matter.
	for (coords_vector::const_iterator it = impurity_coords.begin(); it != impurity_coords.end(); ++it) {
		const double x = std::tr1::get<0>(*it);
		const double y = std::tr1::get<1>(*it);
		const double z = std::tr1::get<2>(*it);
		for (size_t sx=0; sx<dl.sizex; sx++) {
			const double px = dl.get_posx(sx);
			for (size_t sy=0; sy<dl.sizey; sy++) {
				const double py = dl.get_posy(sy);
				const double rx = x-px;
				const double ry = y-py;
				const double r2 = rx*rx + ry*ry + z*z;
				dl.value(pot_values, sx, sy) += alpha*pow(r2, -exponent/2);
			}
		}
	}
}

void CoulombImpurities::init() {
	std::stringstream ss;
	ss << "Coulomb-like impurities, density = " << density << ", exponent = "
		<< exponent << ", strength = " << alpha << ", max displacement = " <<
		maxd;
	description = ss.str();
}
