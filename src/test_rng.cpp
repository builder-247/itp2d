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

// Test that the normally distributed random numbers are really normally distributed.
// This could be done more thoroughly, but here we just check the mean and
// variance of the distribution.
TEST(rng, gaussianity_of_gaussian_rand) {
	const int N = 100000;
	const double tolerance = 0.015;
	RNG rng;
	double* membuf = new double[N];
	for (int i=0; i<N; i++) {
		membuf[i] = rng.gaussian_rand();
	}
	double sum;
	sum = 0;
	for (int i=0; i<N; i++) {
		sum += membuf[i];
	}
	const double mean = sum/N;
	sum = 0;
	for (int i=0; i<N; i++) {
		sum += (mean-membuf[i])*(mean-membuf[i]);
	}
	const double variance = sum/N;
	ASSERT_LT(fabs(mean), tolerance);
	ASSERT_NEAR(variance, 1.0, 100*tolerance);
	if (dump_data) {
		// Save the random numbers for more careful analysis.
		const char filename[] = "data/test_rng_gaussian.h5";
		const hsize_t hN = N;
		H5::H5File datafile(filename, H5F_ACC_TRUNC);
		H5::DataSpace rand_space(1, &hN);
		H5::DataSet rand_data = datafile.createDataSet("rand", H5::PredType::NATIVE_DOUBLE, rand_space);
		rand_data.write(membuf, H5::PredType::NATIVE_DOUBLE);
		datafile.close();
	}
	delete[] membuf;
}

// The same thing for the uniform distribution.
TEST(rng, uniformity_of_uniform_rand) {
	const int N = 100000;
	const double tolerance = 0.015;
	RNG rng;
	double* membuf = new double[N];
	for (int i=0; i<N; i++) {
		membuf[i] = rng.uniform_rand();
	}
	double sum;
	sum = 0;
	for (int i=0; i<N; i++) {
		sum += membuf[i];
	}
	const double mean = sum/N;
	sum = 0;
	for (int i=0; i<N; i++) {
		sum += (mean-membuf[i])*(mean-membuf[i]);
	}
	const double variance = sum/N;
	ASSERT_NEAR(variance, 0.0, 100*tolerance);
	if (dump_data) {
		// Save the random numbers for more careful analysis.
		const char filename[] = "data/test_rng_uniform.h5";
		const hsize_t hN = N;
		H5::H5File datafile(filename, H5F_ACC_TRUNC);
		H5::DataSpace rand_space(1, &hN);
		H5::DataSet rand_data = datafile.createDataSet("rand", H5::PredType::NATIVE_DOUBLE, rand_space);
		rand_data.write(membuf, H5::PredType::NATIVE_DOUBLE);
		datafile.close();
	}
	delete[] membuf;
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
	RNG rng;
	int n = 0;
	for (int i=0; i<N; i++) {
		if (rng.bernoulli_trial(p))
			n++;
	}
	ASSERT_NEAR(p*N, n, tolerance*N);
}

INSTANTIATE_TEST_CASE_P(rng, BernoulliTest, testing::Range(0.1, 0.9, 0.1));
