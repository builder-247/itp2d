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
 * A simple timer using POSIX sys/time.h or mach_absolute_time
 *
 * This is one place that might need work if itp2d is to be ported to non-POSIX
 * platforms (=Windows). However, who runs numerical simulations on Windows?
 */

#ifndef _TIMER_HPP_
#define _TIMER_HPP_

// For OS X we use mach_absolute_time, for POSIX we use clock_gettime
#ifdef __MACH__
#include <mach/mach.h>
#include <mach/mach_time.h>
#else
#include <sys/time.h>
#endif

#ifndef NSEC_PER_SEC
#define NSEC_PER_SEC 1000000000L
#endif

#include "exceptions.hpp"

class Timer {
public:
	#ifdef __MACH__
	typedef uint64_t time;
	#else
	typedef timespec time;
	#endif
	Timer();
	inline void start();
	inline void stop();
	inline void reset();
	inline double get_time();
private:
	time start_time;
	time stop_time;
	unsigned long int elapsed_nsec;
	bool running;
	#ifdef __MACH__
	mach_timebase_info_data_t timebase_info;
	#endif
};

inline void Timer::start() {
	if (running)
		throw GeneralError("Timer::start() called on a timer which was already running");
	#ifdef __MACH__
	start_time = mach_absolute_time();
	#else
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	#endif
	running = true;
}

inline void Timer::stop() {
	if (not running)
		throw GeneralError("Timer::stop() called on a timer which was not running");
	#ifdef __MACH__
	stop_time = mach_absolute_time();
	// convert difference to nanoseconds
	elapsed_nsec += (stop_time-start_time) * timebase_info.numer / timebase_info.denom;
	#else
	clock_gettime(CLOCK_MONOTONIC, &stop_time);
	// Compute difference in nanoseconds
	const long int sec_diff = stop_time.tv_sec - start_time.tv_sec;
	const long int nsec_diff = stop_time.tv_nsec - start_time.tv_nsec;
	elapsed_nsec += (nsec_diff < 0)?
		NSEC_PER_SEC*(sec_diff-1) + (NSEC_PER_SEC + nsec_diff) :
		NSEC_PER_SEC*(sec_diff) + nsec_diff;
	#endif
	running = false;
}

inline void Timer::reset() {
	elapsed_nsec = 0;
	running = false;
}

inline double Timer::get_time() {
	if (running) {
		throw GeneralError("Timer::get_time() called on a timer which was still running");
	}
	return 1e-9*static_cast<double>(elapsed_nsec);
}

#endif // _TIMER_HPP_
