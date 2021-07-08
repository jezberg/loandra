#ifndef MAXPP_CLAUSE_HPP
#define MAXPP_CLAUSE_HPP

#include <vector>
#include <cstdint>

namespace maxPreprocessor {
class ClauseIter;
class ConstClauseIter;

class Clause {
public:
	std::vector<int> lit;
	uint64_t weight;
	uint64_t hash;
	
	void updateHash();
	bool isHard() const;
	void addLiteral(int lit);
	void removeLiteral(int lit);
	Clause (std::vector<int> literals_, uint64_t weight_);
};
}
#endif