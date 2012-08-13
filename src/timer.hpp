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
 * A simple timer using POSIX sys/time.h
 *
 * This is one place that might need work if itp2d is to be ported to non-POSIX
 * platforms (=Windows). However, who runs numerical simulations on Windows?
 */

#ifndef _TIMER_HPP_
#define _TIMER_HPP_

#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
	Timer();
	typedef struct timeval timeval;
	inline void start();	// Starts the timer
	inline double stop();	// Stops and returns elapsed time
private:
	timeval start_time;
	timeval stop_time;
	bool running;
};

inline void Timer::start() {
	gettimeofday(&start_time, NULL);
	running = true;
}

inline double Timer::stop() {
	gettimeofday(&stop_time, NULL);
	if (!running)
		return 0;
	running = false;
	timeval diff;
	timersub(&stop_time, &start_time, &diff);
	return static_cast<double>(diff.tv_sec) + 1e-6*static_cast<double>(diff.tv_usec);
}

#endif // _TIMER_HPP_
