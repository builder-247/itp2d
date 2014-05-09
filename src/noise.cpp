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

// noise parser function

Noise const* parse_noise_description(std::string const& str,
		DataLayout const& dl, Constraint const& constraint, RNG& rng) {
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
		return new NoNoise();
	else if (name == "gaussian" or name == "gaussians" or name == "gaussiannoise") {
		if (params.size() == 3)
			return new GaussianNoise(params[0], params[1], 0.0, params[2], 0.0, dl, constraint, rng);
		else if (params.size() == 5)
			return new GaussianNoise(params[0], params[1], params[2], params[3], params[4], dl, constraint, rng);
		else
			throw InvalidNoiseType("Noise type GaussianNoise takes either 3 or 5 parameters");
	}
	else if (name == "coulomb" or name == "coulombimpurities") {
	// It looks like this "design pattern" has reached its end
		if (params.size() == 4)
			return new CoulombImpurities(params[0], params[1], params[2], params[3], dl, constraint, rng);
		else
			throw InvalidNoiseType("Noise type CoulombImpurities takes 4 parameters");
	}
	else if (name == "hemisphere" or name == "hemispheres" or name == "hemisphereimpurities") {
		if (params.size() == 3)
			return new HemisphereImpurities(params[0], params[1], params[2], dl, constraint, rng);
		else
			throw InvalidNoiseType("Noise type HemisphereImpurities takes 3 parameters");
	}
	else
		throw UnknownNoiseType(str);
}

// GaussianNoise

GaussianNoise::GaussianNoise(double d, double amp, double amp_stdev, double w, double w_stdev,
		DataLayout const& dl, Constraint const& constr, RNG& rng) :
		datalayout(dl), constraint(constr),
		density(d), amplitude_mean(amp), amplitude_stdev(amp_stdev), width_mean(w), width_stdev(w_stdev) {
	init(rng);
}

void GaussianNoise::init(RNG& rng) {
	// Write description
	std::stringstream ss;
	ss << "gaussian spikes, density = " << density << ", amplitude ~ N(" << amplitude_mean
		<< "," << amplitude_stdev << "^2), width ~ N(" << width_mean << "," << width_stdev << "^2)";
	description = ss.str();
	// Generate noise
	spikes.clear();
	// The number of impurities on the calculation plane is Poisson distributed around the mean value
	const double lambda = density*datalayout.lenx*datalayout.leny;
	const unsigned int N = rng.poisson_rand(lambda);
	spikes.reserve(N);
	// First, randomly distribute the centers, amplitudes and widths of the spikes
	for (unsigned int n=0; n<N; n++) {
		// First randomize the position
		const double x = (rng.uniform_rand()-0.5)*datalayout.lenx;
		const double y = (rng.uniform_rand()-0.5)*datalayout.leny;
		if (not constraint.check(x, y)) {
			continue;
		}
		// Determine amplitude and width
		const double A = amplitude_mean + amplitude_stdev*rng.gaussian_rand();
		const double w = width_mean + width_stdev*rng.gaussian_rand();
		if (w < 0)
			continue;
		spikes.push_back(spike(x, y, A, w));
	}
}

void GaussianNoise::add_noise(DataLayout const& dl, double* pot_values) const {
	// Loop through all spikes and add their effect into pot_values. This is
	// not a very efficient way to do this, but this only needs to be done once
	// so it doesn't really matter.
	for (std::vector<spike>::const_iterator it = spikes.begin(); it != spikes.end(); ++it) {
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

void GaussianNoise::write_realization_data(std::vector<double>& vec) const {
	vec.clear();
	vec.reserve(4*spikes.size());
	for (std::vector<spike>::const_iterator it = spikes.begin(); it != spikes.end(); ++it) {
		vec.push_back(std::tr1::get<0>(*it));
		vec.push_back(std::tr1::get<1>(*it));
		vec.push_back(std::tr1::get<2>(*it));
		vec.push_back(std::tr1::get<3>(*it));
	}
}

// CoulombImpurities

CoulombImpurities::CoulombImpurities(double d, double e, double a, double _maxd, DataLayout const& dl, Constraint const& constr, RNG& rng) :
		datalayout(dl), constraint(constr),
		density(d), exponent(e), alpha(a), maxd(_maxd) {
	init(rng);
}

void CoulombImpurities::init(RNG& rng) {
	// Write description
	std::stringstream ss;
	ss << "Coulomb-like impurities, density = " << density << ", exponent = "
		<< exponent << ", strength = " << alpha << ", max displacement = " <<
		maxd;
	description = ss.str();
	// Generate noise
	impurities.clear();
	// The number of impurities inside the box is Poisson distributed around the mean value
	const double lambda = density*datalayout.lenx*datalayout.leny*2*maxd;
	const unsigned int N = rng.poisson_rand(lambda);
	impurities.reserve(N);
	// First, randomly distribute the positions of the impurities
	for (unsigned int n=0; n<N; n++) {
		const double x = (rng.uniform_rand()-0.5)*datalayout.lenx;
		const double y = (rng.uniform_rand()-0.5)*datalayout.leny;
		if (not constraint.check(x, y)) {
			continue;
		}
		const double z = (rng.uniform_rand()-0.5)*2*maxd;
		impurities.push_back(impurity(x, y, z, alpha));
	}
}

void CoulombImpurities::add_noise(DataLayout const& dl, double* pot_values) const {
	for (std::vector<impurity>::const_iterator it = impurities.begin(); it != impurities.end(); ++it) {
		const double x = std::tr1::get<0>(*it);
		const double y = std::tr1::get<1>(*it);
		const double z = std::tr1::get<2>(*it);
		const double a = std::tr1::get<3>(*it);
		for (size_t sx=0; sx<dl.sizex; sx++) {
			const double px = dl.get_posx(sx);
			for (size_t sy=0; sy<dl.sizey; sy++) {
				const double py = dl.get_posy(sy);
				const double rx = x-px;
				const double ry = y-py;
				const double r2 = rx*rx + ry*ry + z*z;
				dl.value(pot_values, sx, sy) += a*pow(r2, -exponent/2);
			}
		}
	}
}

void CoulombImpurities::write_realization_data(std::vector<double>& vec) const {
	vec.clear();
	vec.reserve(4*impurities.size());
	for (std::vector<impurity>::const_iterator it = impurities.begin(); it != impurities.end(); ++it) {
		vec.push_back(std::tr1::get<0>(*it));
		vec.push_back(std::tr1::get<1>(*it));
		vec.push_back(std::tr1::get<2>(*it));
		vec.push_back(std::tr1::get<3>(*it));
	}
}

// HemisphereImpurities

HemisphereImpurities::HemisphereImpurities(double d, double A, double r, DataLayout const& dl, Constraint const& constr, RNG& rng) :
		datalayout(dl), constraint(constr),
		density(d), amplitude(A), radius(r) {
	init(rng);
}

void HemisphereImpurities::init(RNG& rng) {
	// Write description
	std::stringstream ss;
	ss << "Hemisphere impurities, density = " << density << ", amplitude = "
		<< amplitude << ", radius = " << radius;
	description = ss.str();
	// Generate noise
	impurities.clear();
	// The number of impurities inside the box is Poisson distributed around the mean value
	const double lambda = density*datalayout.lenx*datalayout.leny;
	const unsigned int N = rng.poisson_rand(lambda);
	impurities.reserve(N);
	// First, randomly distribute the positions of the impurities
	for (unsigned int n=0; n<N; n++) {
		const double x = (rng.uniform_rand()-0.5)*datalayout.lenx;
		const double y = (rng.uniform_rand()-0.5)*datalayout.leny;
		if (not constraint.check(x, y)) {
			continue;
		}
		impurities.push_back(impurity(x, y, amplitude, radius));
	}
}

void HemisphereImpurities::add_noise(DataLayout const& dl, double* pot_values) const {
	for (std::vector<impurity>::const_iterator it = impurities.begin(); it != impurities.end(); ++it) {
		const double x = std::tr1::get<0>(*it);
		const double y = std::tr1::get<1>(*it);
		const double A = std::tr1::get<2>(*it);
		const double r = std::tr1::get<3>(*it);
		const double r2 = r*r;
		for (size_t sx=0; sx<dl.sizex; sx++) {
			const double px = dl.get_posx(sx);
			for (size_t sy=0; sy<dl.sizey; sy++) {
				const double py = dl.get_posy(sy);
				const double dx = x-px;
				const double dy = y-py;
				const double d2 = dx*dx + dy*dy;
				if (d2 < r2)
					dl.value(pot_values, sx, sy) += A*sqrt(1-d2/r2);
			}
		}
	}
}

void HemisphereImpurities::write_realization_data(std::vector<double>& vec) const {
	vec.clear();
	vec.reserve(4*impurities.size());
	for (std::vector<impurity>::const_iterator it = impurities.begin(); it != impurities.end(); ++it) {
		vec.push_back(std::tr1::get<0>(*it));
		vec.push_back(std::tr1::get<1>(*it));
		vec.push_back(std::tr1::get<2>(*it));
		vec.push_back(std::tr1::get<3>(*it));
	}
}
