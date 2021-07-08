#ifndef MAXPP_PROBLEMINSTANCE_HPP
#define MAXPP_PROBLEMINSTANCE_HPP

#include <vector>
#include <cstdint>

#include "clause.hpp"
#include "global.hpp"
#include "touchedlist.hpp"

namespace maxPreprocessor {
class ProblemInstance {
public:
	
	// All clauses including removed. Clauses id is its index in this list and it never changes.
	std::vector<Clause> clauses;
	std::vector<std::vector<int> > litClauses;
	// Is clause i removed
	std::vector<int> removedClauses;
	
	// Is variable i a label. 0 if not, VAR_FALSE if its negation is soft, VAR_TRUE otherwise
	std::vector<int> isLabel;
	
	int vars;
	int excessVar;
	
	TouchedList tl;
	
	ProblemInstance(const std::vector<std::vector<int> >& clauses_, const std::vector<uint64_t>& weights_, uint64_t topWeight);
	
	bool isSimpleSoftClause(int clause) const;
	bool isClauseRemoved(int clauseId) const;
	bool isHard(int clauseId) const;
	bool isVarRemoved(int var) const;
	uint64_t labelWeight(int lbl) const;
	
	void populateLitClauses(int clauseId);
	
	void removeClauseFromLitClause(int clauseId, int lit);
	void removeClauseFromLitClauses(int clauseId);
	
	void removeClause(int clauseId);
	
	// Literals in clause should be in sorted order
	void addClause(const std::vector<int>& clause, uint64_t weight = HARDWEIGHT);
	
	// Returns the index of added variable
	int addVar();
	
	void addLiteralToClause(int lit, int clause, bool touch = true);
	
	// The clause must be hard
	void removeLiteralFromClause(int lit, int clause, bool touch = true);
	
	bool canSubsume2(uint64_t h1, uint64_t h2);
	bool canSubsume1(int clause1, int clause2);
	bool canSubsume(int clause1, int clause2);
	
	int getExcessVar();
	uint64_t getWeightSum();
};
}
#endif