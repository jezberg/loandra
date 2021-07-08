#ifndef MAXPP_PREPROCESSEDINSTANCE_HPP
#define MAXPP_PREPROCESSEDINSTANCE_HPP

#include <vector>
#include <iostream>
#include <cstdint>

namespace maxPreprocessor {
class PreprocessedInstance {
public:
	std::vector<std::vector<int> > clauses;
	std::vector<uint64_t> weights;
	std::vector<std::pair<int, uint64_t> > labels;
};
}
#endif