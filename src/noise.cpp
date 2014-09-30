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

using namespace std;

// GaussianImpurities

void GaussianImpurities::add_noise(double sx, double sy, vector<double> const& params, DataLayout const& dl, double* pot_values) const {
	if (params.size() != num_params)
		throw GeneralError("GaussianImpurities constructor with incorrect length for parameter vector. This should never happen.");
	const double A = params[0];
	const double w = params[1];
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

// SpatialImpurities

SpatialImpurities::SpatialImpurities(ImpurityType const& _type, ImpurityDistribution const& _distribution) :
		type(_type), distribution(_distribution) {
	list<coordinate_pair> const& coords = distribution.get_coordinates();
	for (list<coordinate_pair>::const_iterator it = coords.begin(); it != coords.end(); ++it) {
		double const& x = it->first;
		double const& y = it->second;
		realization_data.push_back(x);
		realization_data.push_back(y);
		type.new_realization(realization_data);
	}
}

void SpatialImpurities::add_noise(DataLayout const& dl, double* pot_values) const {
	if (realization_data.size() % (2+type.get_num_params()) != 0) {
		throw GeneralError("SpatialImpurities::add_noise() with incorrect length for realization_data vector. This should never happen.");
	}
	list<double>::const_iterator it = realization_data.begin();
	vector<double> params;
	do {
		const double x = *it; ++it;
		const double y = *it; ++it;
		for (size_t i=0; i<type.get_num_params(); i++) {
			params.push_back(*it); it++;
		}
		type.add_noise(x, y, params, dl, pot_values);
		params.clear();
	} while (it != realization_data.end());
}

void SpatialImpurities::write_realization_data(vector<double>& vec) const {
	vec.clear();
	vec.reserve(realization_data.size());
	copy(realization_data.begin(), realization_data.end(), back_inserter(vec));
}

// noise parser function

ImpurityType const* parse_impurity_type_description(string const& type_str, DataLayout const& dl, RNG& rng) {
	name_parameters_pair p;
	p = parse_parameter_string(type_str);
	string const& name = p.first;
	vector<double> const& params = p.second;
	// Simply delegate to the individual constructors based on name
	if (name == "gaussian" or name == "gaussians" or name == "gaussiannoise") {
		if (params.size() == 2)
			return new GaussianImpurities(params[0], 0.0, params[1], 0.0, rng);
		else if (params.size() == 4)
			return new GaussianImpurities(params[0], params[1], params[2], params[3], rng);
		else
			throw InvalidImpurityType("Impurity type GaussianImpurities takes either 2 or 4 parameters");
	}
	else
		throw UnknownImpurityType(type_str);
}

ImpurityDistribution const* parse_impurity_distribution_description(string const& distribution_str, DataLayout const& dl, Constraint const& constraint, RNG& rng) {
	name_parameters_pair p;
	p = parse_parameter_string(distribution_str);
	string const& name = p.first;
	vector<double> const& params = p.second;
	// Simply delegate to the individual constructors based on name
	if (name == "uniform") {
		if (params.size() == 1)
			return new UniformImpurities(params[0], dl, constraint, rng);
		else
			throw InvalidImpurityDistributionType("Impurity disribution type UniformImpurities takes only one parameter");
	}
	else
		throw UnknownImpurityDisributionType(distribution_str);
}

/*

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

// HemisphereImpurities

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

*/
