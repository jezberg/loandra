#ifndef MAXPP_TRACE_HPP
#define MAXPP_TRACE_HPP

#include <vector>
#include <iostream>
#include <cstdint>

namespace maxPreprocessor {
class Trace {
public:
	uint64_t removedWeight = 0;
	std::vector<int> operations;
	std::vector<std::vector<int> > data;
	
	void setVar(int var, bool value);
	void BVE(int var, const std::vector<std::vector<int> >& nClauses);
	void BCE(int lit, const std::vector<int>& clause);
	void BCR(int lbl1, int lbl2, const std::vector<std::vector<int> >& nClauses);
	void setEqual(int lit1, int lit2);
	void LS(int lbl, int lit, const std::vector<int>& clause);
	void labelEliminate(int lbl1, int lbl2, int tautli);
	void removeWeight(uint64_t weight);
	std::pair<std::vector<int>, uint64_t> getSolution(const std::vector<int>& trueLits, uint64_t weight, int vars, int originalVars);
	void printSolution(std::ostream& output, const std::vector<int>& trueLits, uint64_t weight, int vars, int originalVars);
};
}
#endif