#include <chrono>

#include "timer.hpp"

using namespace std;
namespace maxPreprocessor {

Timer::Timer() {
	timing = false;
	elapsedTime = chrono::duration<double>(std::chrono::duration_values<double>::zero());
}

void Timer::start() {
	if (timing) return;
	timing = true;
	startTime = chrono::steady_clock::now();
}

void Timer::stop() {
	if (!timing) return;
	timing = false;
	chrono::time_point<std::chrono::steady_clock> endTime = chrono::steady_clock::now();
	elapsedTime += (endTime - startTime);
}

chrono::duration<double> Timer::getTime() {
	if (timing) {
		stop();
		chrono::duration<double> ret = elapsedTime;
		start();
		return ret;
	}
	else {
		return elapsedTime;
	}
}

}