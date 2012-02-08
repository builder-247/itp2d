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

/* Interface class representing "noise" which can be added to some array of
 * values, such as a potential. */
class Noise {
	public:
		virtual ~Noise() {};
		virtual void add_noise(DataLayout const& dl, double* pot_values, RNG& rng) const = 0;
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
		NoNoise() {
			init();
		}
		NoNoise(std::vector<double> params) {
			if (not params.empty())
				throw InvalidNoiseType("Noise type NoNoise does not take parameters");
			init();
		}
		void add_noise(__attribute__((unused)) DataLayout const& dl, __attribute__((unused))double* pot_values, __attribute__((unused))RNG& rng) const {}
	private:
		void init() {
			description = "none";
		}
};

/* Randomly distributed gaussian spikes with normally distributed width
 * and normally distributed amplitude. */
class GaussianSpikes : public Noise {
	public:
		static const double default_relative_amplitude_stdev;
		static const double default_relative_width_stdev;
		GaussianSpikes(double _density, double _amplitude_mean, double _width_mean) :
				density(_density), amplitude_mean(_amplitude_mean), width_mean(_width_mean) {
				amplitude_stdev = default_relative_amplitude_stdev * amplitude_mean;
				width_stdev = default_relative_width_stdev * width_mean;
			init();
		}
		GaussianSpikes(double _density, double _amplitude_mean, double _amplitude_stdev, double _width_mean, double _width_stdev) :
				density(_density), amplitude_mean(_amplitude_mean),
				amplitude_stdev(_amplitude_stdev), width_mean(_width_mean),
				width_stdev(_width_stdev) {
			init();
		}
		GaussianSpikes(std::vector<double> params) {
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
				throw InvalidNoiseType("Noise type GaussianSpikes takes either 3 or 5 parameters");
			init();
		}
		void add_noise(DataLayout const& dl, double* pot_values, RNG& rng) const;
	private:
		double density;
		double amplitude_mean;
		double amplitude_stdev;
		double width_mean;
		double width_stdev;
		void init() {
			std::stringstream ss;
			ss << "gaussian spikes, density = " << density << ", amplitude ~ N(" << amplitude_mean
				<< "," << amplitude_stdev << "^2), width ~ N(" << width_mean << "," << width_stdev << "^2)";
			description = ss.str();
		}
};

#endif // _NOISE_HPP_
