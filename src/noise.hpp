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

#ifndef _NOISE_HPP_
#define _NOISE_HPP_

#include <vector>
#include <tr1/tuple>

#include "itp2d_common.hpp"
#include "datalayout.hpp"
#include "rng.hpp"
#include "parser.hpp"
#include "constraint.hpp"

/* Interface class representing "noise" which can be added to some array of
 * values, such as a potential. */
class Noise {
	public:
		virtual ~Noise() {};
		virtual void add_noise(DataLayout const& dl, double* pot_values, RNG& rng, Constraint const& constraint) const = 0;
		virtual std::string const& get_description() const { return description; }
	protected:
		std::string description;
};

// A parser function for returning a Noise instance from a user provided
// desciption string

Noise const* parse_noise_description(std::string const& str);

// Individual noise types
// A noise type must implement
// 	* The function add_noise for adding noise to an array of potential values
// 	* A constructor constructing the instance from a vector of doubles. This is
// 	  used for constructing noise types from user-provided strings

/* The simple case of no noise at all */
class NoNoise : public Noise {
	public:
		NoNoise() { init(); }
		NoNoise(std::vector<double> params);
		void add_noise(__attribute__((unused)) DataLayout const& dl, __attribute__((unused))double* pot_values, __attribute__((unused))RNG& rng, __attribute__((unused)) Constraint const& constraint) const {}
	private:
		void init() { description = "none"; }
};

/* Randomly distributed gaussian spikes with normally distributed width
 * and normally distributed amplitude. */
class GaussianNoise : public Noise {
	public:
		static const double default_relative_amplitude_stdev;
		static const double default_relative_width_stdev;
		GaussianNoise(double _density, double _amplitude_mean, double _width_mean) :
				density(_density), amplitude_mean(_amplitude_mean), width_mean(_width_mean) {
				amplitude_stdev = default_relative_amplitude_stdev * amplitude_mean;
				width_stdev = default_relative_width_stdev * width_mean;
			init();
		}
		GaussianNoise(double _density, double _amplitude_mean, double _amplitude_stdev, double _width_mean, double _width_stdev) :
				density(_density), amplitude_mean(_amplitude_mean),
				amplitude_stdev(_amplitude_stdev), width_mean(_width_mean),
				width_stdev(_width_stdev) {
			init();
		}
		GaussianNoise(std::vector<double> params);
		void add_noise(DataLayout const& dl, double* pot_values, RNG& rng, Constraint const& constraint) const;
	private:
		double density;
		double amplitude_mean;
		double amplitude_stdev;
		double width_mean;
		double width_stdev;
		void init();
};

/* Coulomb-like impurities in 3D space with a tunable exponent. */
class CoulombImpurities : public Noise {
	public:
		static const double default_alpha;
		static const double default_exponent;
		static const double default_maxd;
		CoulombImpurities(double _density, double _exponent = default_exponent,
				double _alpha = default_alpha, double _maxd = default_maxd) :
				density(_density), exponent(_exponent), alpha(_alpha), maxd(_maxd) {
			init();
		}
		CoulombImpurities(std::vector<double> params);
		void add_noise(DataLayout const& dl, double* pot_values, RNG& rng, Constraint const& constraint) const;
	private:
		double density; // density of impurities in three dimensions
		double exponent;
		double alpha;	// the potential of a single impurity is alpha/r^exp
		double maxd;	// the impurities are located furthest at maxd distance units from the calculation plane
		void init();
};

#endif // _NOISE_HPP_
