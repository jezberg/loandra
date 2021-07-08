#ifndef MAXPP_TOUCHEDLIST_HPP
#define MAXPP_TOUCHEDLIST_HPP
#include <vector>
#include <deque>
#include <cstdint>
#include <map>

namespace maxPreprocessor {
class ProblemInstance;

class TouchedList {
private:
	std::deque<std::pair<uint64_t, int> > touchedLiterals;
	std::deque<std::pair<uint64_t, int> > modLiterals;
	std::deque<std::pair<uint64_t, int> > touchedClauses;
	std::deque<std::pair<uint64_t, int> > modClauses;
	std::map<std::string, uint64_t> techniqueItr;
	void getTouchedLiteralsCh(uint64_t eitr, std::vector<int>& ret);
	void getModLiteralsCh(uint64_t eitr, std::vector<int>& ret);
	void getTouchedVarsCh(uint64_t eitr, std::vector<int>& ret);
	void getModVarsCh(uint64_t eitr, std::vector<int>& ret);
	std::vector<int> getI;
	uint64_t itr;
	uint64_t frontItr;
	int getItr;
	int vars;
	const ProblemInstance& pi;
public:
	TouchedList(const ProblemInstance& pi_);
	void init(int vars_);
	
	void shrink(bool log);
	
	void addVar();
	void modClause(int c);
	void modLiteral(int l);
	void touchClause(int c);
	void touchLiteral(int l);
	
	void setItr(std::string technique);
	
	std::vector<int> getModClauses(std::string technique);
	std::vector<int> getModLiterals(std::string technique);
	std::vector<int> getModVariables(std::string technique);
	std::vector<int> getTouchedLiterals(std::string technique);
	std::vector<int> getTouchedVariables(std::string technique);
	std::vector<int> getBinaryLiterals(std::string technique);
};
}

#endif