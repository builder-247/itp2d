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
 * Unit tests for the random number generator.
 */

#include "test_rng.hpp"

// Helper function for computing the mean and variance of a sequence of numbers
std::pair<double,double> mean_and_variance(std::vector<double> const& vec) {
	const size_t N = vec.size();
	double sum;
	sum = 0;
	for (size_t i=0; i<N; i++) {
		sum += vec[i];
	}
	const double mean = sum/static_cast<double>(N);
	sum = 0;
	for (size_t i=0; i<N; i++) {
		sum += (mean-vec[i])*(mean-vec[i]);
	}
	const double variance = sum/static_cast<double>(N);
	return std::make_pair(mean, variance);
}

void write_sample(std::vector<double> const& vec, std::string filename) {
	const hsize_t hN = vec.size();
	H5::H5File datafile(filename.c_str(), H5F_ACC_TRUNC);
	H5::DataSpace rand_space(1, &hN);
	H5::DataSet rand_data = datafile.createDataSet("rand", H5::PredType::NATIVE_DOUBLE, rand_space);
	rand_data.write(&vec.front(), H5::PredType::NATIVE_DOUBLE);
	datafile.close();
}

// Test that the normally distributed random numbers are really normally distributed.
// This could be done more thoroughly, but here we just check the mean and
// variance of the distribution.
TEST(rng, gaussianity_of_gaussian_rand) {
	const size_t N = 100000;
	const double tolerance = 0.01;
	RNG rng(RNG::produce_random_seed());
	std::vector<double> sample(N);
	for (size_t i=0; i<N; i++) {
		sample[i] = rng.gaussian_rand();
	}
	std::pair<double,double> p = mean_and_variance(sample);
	const double mean = p.first;
	const double variance = p.second;
	EXPECT_NEAR(mean, 0.0, tolerance);
	EXPECT_NEAR(variance, 1.0, 2*tolerance);
	if (dump_data) {
		// Save the random numbers for more careful analysis.
		write_sample(sample, "data/test_rng_gaussian.h5");
	}
}

// The same thing for the uniform distribution.
TEST(rng, uniformity_of_uniform_rand) {
	const size_t N = 100000;
	const double tolerance = 0.01;
	RNG rng(RNG::produce_random_seed());
	std::vector<double> sample(N);
	for (size_t i=0; i<N; i++) {
		sample[i] = rng.uniform_rand();
	}
	std::pair<double,double> p = mean_and_variance(sample);
	const double mean = p.first;
	const double variance = p.second;
	EXPECT_NEAR(mean, 0.5, tolerance);
	EXPECT_NEAR(variance, 1.0/12, 2*tolerance);
	if (dump_data) {
		write_sample(sample, "data/test_rng_uniform.h5");
	}
}

// The same thing for the Poisson distribution.
TEST(rng, poissonity_of_poisson_rand) {
	const double lambda = 1.0;
	const size_t N = 100000;
	const double tolerance = 0.01;
	RNG rng(RNG::produce_random_seed());
	std::vector<double> sample(N);
	for (size_t i=0; i<N; i++) {
		sample[i] = rng.poisson_rand(lambda);
	}
	std::pair<double,double> p = mean_and_variance(sample);
	const double mean = p.first;
	const double variance = p.second;
	EXPECT_NEAR(mean, lambda, tolerance);
	EXPECT_NEAR(variance, lambda, 2*tolerance);
	if (dump_data) {
		write_sample(sample, "data/test_rng_uniform.h5");
	}
}

class BernoulliTest : public ::testing::TestWithParam<double> {
	public:
		BernoulliTest() : p(GetParam()) {}
		const double p;
};

// Test the mean of the bernoulli distribution for some range of parameters.
TEST_P(BernoulliTest, is_bernoulli) {
	const int N = 100000;
	const double tolerance = 0.01;
	RNG rng(RNG::produce_random_seed());
	int n = 0;
	for (int i=0; i<N; i++) {
		if (rng.bernoulli_trial(p))
			n++;
	}
	EXPECT_NEAR(p*N, n, tolerance*N);
}

INSTANTIATE_TEST_CASE_P(rng, BernoulliTest, testing::Range(0.1, 0.9, 0.1));
