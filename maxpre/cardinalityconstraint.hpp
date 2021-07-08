#ifndef MAXPP_CARDINALITYCONSTRAINT_HPP
#define MAXPP_CARDINALITYCONSTRAINT_HPP

#include <cstdint>
#include <vector>

#include "global.hpp"

namespace maxPreprocessor {
class CardinalityConstraint {
private:
	struct Formula {
		std::vector<std::vector<int> > clauses;
		int freeVar;
		size_t size() const;
		Formula (int freeVar_);
	};
	std::vector<int> SMerge(const std::vector<int>& A, const std::vector<int>& B, Formula& F, bool d1, bool d2);
	std::vector<int> HMerge(const std::vector<int>& A, const std::vector<int>& B, Formula& F, bool d1, bool d2);
	std::vector<int> HSort(const std::vector<int>& A, int l, int r, Formula& F, bool d1, bool d2);
	std::vector<int> Card(const std::vector<int>& A, int l, int r, int k, Formula& F, bool d1, bool d2);
	std::pair<Formula, int> getLessNetwork(int freeVar, std::vector<int> input);
	std::pair<Formula, int> getGreaterNetwork(int freeVar, std::vector<int> input);
	std::pair<Formula, int> getBestLessNetwork(int freeVar, std::vector<int> input);
	std::pair<Formula, int> getEqualNetwork(int freeVar, std::vector<int> input);
	std::pair<Formula, int> getNequalNetwork(int freeVar, std::vector<int> input);
public:
	static const int CONSTRAINT_EQUAL = 1;
	static const int CONSTRAINT_NEQUAL = 2;
	static const int CONSTRAINT_LESS_EQ = 3;
	static const int CONSTRAINT_ALL_LESS = 4;
	static const int CONSTRAINT_ALL_GREATER = 5;
	static const int CONSTRAINT_ALL_BOTH = 6;
	
	static const int CODING_NETWORK = 1;
	static const int CODING_DEFAULT = CODING_NETWORK;
	
	const int type;
	const int coding;
	const std::vector<int> inputLiterals;
	const std::vector<int> outputLiterals;
	const int K;
	const std::vector<int> addLiterals;
	const bool iff;
	
	std::vector<std::vector<int> > getClauses(int freeVar);
	CardinalityConstraint(int type_, int coding_, const std::vector<int>& inputLiterals_, const std::vector<int>& outputLiterals_, int K_, const std::vector<int>& addLiterals_, bool iff_);
};
}
#endif