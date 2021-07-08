#ifndef MAXPP_GLOBAL_HPP
#define MAXPP_GLOBAL_HPP

#include <cstdint>

namespace maxPreprocessor {
const uint64_t HARDWEIGHT = (uint64_t)1 << (uint64_t)63;

const int VAR_UNDEFINED = 0;
const int VAR_TRUE = 1;
const int VAR_FALSE = 2;

const int INPUT_FORMAT_WPMS = 1;
const int INPUT_FORMAT_MS = 2;
const int INPUT_FORMAT_SAT = 3;

inline int litValue(int x) {
	return !(x&1);
}
inline int litVariable(int x) {
	return (x>>1);
}
inline int litNegation(int x) {
	return x^1;
}

inline int posLit(int x) {
	return x<<1;
}
inline int negLit(int x) {
	return (x<<1)^1;
}
}
#endif