#ifndef MAXPP_AMSLEX_HPP
#define MAXPP_AMSLEX_HPP

#include <vector>
#include <map>

#include "probleminstance.hpp"
#include "timer.hpp"

namespace maxPreprocessor{
class AMSLEX{
private:
	std::vector<int> ALU;
	std::vector<int> ALI;
	int ALIt;
	std::vector<int> ALF;
	const unsigned BSconst = 18; // magic constant
	
	struct vecP {
		unsigned B, E, O;
		unsigned size() const {
			return E - B;
		}
	};
	
	int* data;
	unsigned dataSize;
	const ProblemInstance& pi;
	bool isPrefix(const std::vector<int>& a, const std::vector<int>& b);
	bool isPrefix(vecP a, vecP b);
	void assumeSize(unsigned size);
	bool CSO1(const std::vector<vecP>& D, unsigned b, unsigned e, const vecP S, 
					   unsigned j, unsigned d, const std::vector<int>& d0Index);
	bool CSO2(const std::vector<int>& D, unsigned b, unsigned e, const std::vector<int>& S, unsigned j, unsigned d);
	std::vector<int> amsLexSEPerm(const std::vector<int>& clauses);
	std::vector<int> amsLexSENonPerm(const std::vector<int>& clauses);
	std::pair<std::vector<int>, std::vector<int> > amsLexSSRPerm(const std::vector<int>& c1, const std::vector<int>& c2, int var);
public:
	AMSLEX(const ProblemInstance& pi_);
	~AMSLEX();
	
	std::vector<int> amsLexSE(const std::vector<int>& clauses);
	std::pair<std::vector<int>, std::vector<int> > amsLexSSR(const std::vector<int>& c1, const std::vector<int>& c2, int var);
};
}
#endif