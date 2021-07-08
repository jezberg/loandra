#ifndef MAXPP_OUTPUTREADER_HPP
#define MAXPP_OUTPUTREADER_HPP

#include <iostream>
#include <vector>
#include <cstdint>

namespace maxPreprocessor {
class OutputReader {
public:
	int status;
	std::vector<int> trueLits;
	uint64_t ansW;
	
	int readSolution(std::istream& input);
};
}
#endif