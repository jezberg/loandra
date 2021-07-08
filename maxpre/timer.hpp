#ifndef MAXPP_TIMER_HPP
#define MAXPP_TIMER_HPP

#include <chrono>

namespace maxPreprocessor {
class Timer {
private:
	bool timing;
	std::chrono::duration<double> elapsedTime;
	std::chrono::time_point<std::chrono::steady_clock> startTime;
public:
	Timer();
	void start();
	void stop();
	std::chrono::duration<double> getTime();
};
}
#endif