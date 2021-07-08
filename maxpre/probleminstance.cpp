#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <random>

#include "probleminstance.hpp"
#include "utility.hpp"
#include "global.hpp"
#include "touchedlist.hpp"

using namespace std;
namespace maxPreprocessor{

bool ProblemInstance::isClauseRemoved(int clause) const {
	return removedClauses[clause];
}

bool ProblemInstance::isSimpleSoftClause(int clause) const {
	if (!clauses[clause].isHard()) return false;
	int label = -1;
	for (int l : clauses[clause].lit) {
		if (isLabel[litVariable(l)]) {
			if (label != -1) return false;
			if (litClauses[l].size() != 1) return false;
			label = l;
		}
	}
	if (label == -1) return false;
	return true;
}

bool ProblemInstance::isVarRemoved(int var) const {
	if (litClauses[negLit(var)].size() + litClauses[posLit(var)].size() > 0) return false;
	return true;
}

uint64_t ProblemInstance::labelWeight(int lbl) const {
	if (isLabel[lbl] == VAR_TRUE) {
		assert(litClauses[posLit(lbl)].size() == 1);
		return clauses[litClauses[posLit(lbl)][0]].weight;
	}
	else if(isLabel[lbl] == VAR_FALSE) {
		assert(litClauses[negLit(lbl)].size() == 1);
		return clauses[litClauses[negLit(lbl)][0]].weight;
	}
	else {
		assert(0);
	}
}

ProblemInstance::ProblemInstance(const vector<vector<int> >& clauses_, const vector<uint64_t>& weights_, uint64_t topWeight) : tl(*this) {
	assert(clauses_.size() == weights_.size());
	excessVar = 0;
	
	int maxVar = 0;
	
	clauses.reserve(clauses_.size());
	for (unsigned i = 0; i < clauses_.size(); i++) {
		clauses.emplace_back(Clause(clauses_[i], weights_[i]));
	}
	
	for (unsigned i = 0; i < clauses.size(); i++) {
		for (int lit : clauses[i].lit) {
			maxVar = max(maxVar, abs(lit) - 1);
		}
	}
	
	for (unsigned i = 0; i < clauses.size(); i++) {
		for (int& lit : clauses[i].lit) {
			lit=litFromDimacs(lit);
		}
		
		if (clauses[i].weight >= topWeight) clauses[i].weight = HARDWEIGHT;
		
		sort(clauses[i].lit.begin(), clauses[i].lit.end());
		clauses[i].lit.erase(unique(clauses[i].lit.begin(), clauses[i].lit.end()), clauses[i].lit.end());
		clauses[i].updateHash();
	}
	vars = maxVar + 1;
	
	isLabel.resize(vars);
	litClauses.resize(vars*2);
	removedClauses.resize(clauses.size());
	
	tl.init(vars);
	
	for (unsigned i = 0; i < clauses.size(); i++) {
		populateLitClauses(i);
	}
}

void ProblemInstance::populateLitClauses(int clause) {
	for (int lit : clauses[clause].lit) {
		litClauses[lit].push_back(clause);
	}
}

void ProblemInstance::removeClauseFromLitClause(int clause, int lit) {
	for (unsigned i = 0; i < litClauses[lit].size(); i++) {
		if (litClauses[lit][i] == clause) {
			litClauses[lit][i] = litClauses[lit][litClauses[lit].size() - 1];
			litClauses[lit].pop_back();
			break;
		}
	}
}

void ProblemInstance::removeClauseFromLitClauses(int clause) {
	for (int lit : clauses[clause].lit) {
		removeClauseFromLitClause(clause, lit);
	}
}

void ProblemInstance::removeClause(int clause) {
	assert(removedClauses[clause] == false);
	
	tl.touchClause(clause);
	
	removedClauses[clause] = true;
	removeClauseFromLitClauses(clause);
}

// Literals in clause should be in sorted order
void ProblemInstance::addClause(const vector<int>& clause, uint64_t weight) {
	clauses.push_back(Clause(clause, weight));
	int cId = clauses.size() - 1;
	populateLitClauses(cId);
	
	removedClauses.push_back(0);
	
	tl.modClause(cId);
}

// Returns the index of added variable
int ProblemInstance::addVar() {
	litClauses.push_back(vector<int>());
	litClauses.push_back(vector<int>());
	isLabel.push_back(0);
	tl.addVar();
	return vars++;
}

void ProblemInstance::addLiteralToClause(int lit, int clause, bool touch) {
	// Dont add if it already exists
	for (int l : clauses[clause].lit) {
		if (l == lit) return;
		// Tautology -> runtime error
		assert(l != litNegation(lit));
	}
	if (touch) tl.modClause(clause);
	clauses[clause].addLiteral(lit);
	litClauses[lit].push_back(clause);
}

// The clause must be hard
void ProblemInstance::removeLiteralFromClause(int lit, int clause, bool touch) {
	if (touch) {
		tl.modClause(clause);
		tl.touchLiteral(lit);
	}
	clauses[clause].removeLiteral(lit);
	removeClauseFromLitClause(clause, lit);
}

bool ProblemInstance::canSubsume2(uint64_t h1, uint64_t h2) {
	uint64_t un = h1 | h2;
	return un == h1 || un == h2;
}

bool ProblemInstance::canSubsume1(int clause1, int clause2) {
	uint64_t un = clauses[clause1].hash | clauses[clause2].hash;
	return un == clauses[clause1].hash || un == clauses[clause2].hash;
}

bool ProblemInstance::canSubsume(int clause1, int clause2) {
	if ((clauses[clause1].hash | clauses[clause2].hash) == clauses[clause2].hash) return true;
	return false;
}

int ProblemInstance::getExcessVar() {
	if (excessVar == 0) excessVar = addVar();
	return excessVar;
}

uint64_t ProblemInstance::getWeightSum() {
	uint64_t weightSum = 0;
	for (int c = 0; c < (int)clauses.size(); c++) {
		if (!clauses[c].isHard() && !isClauseRemoved(c)) {
			weightSum += clauses[c].weight;
		}
	}
	return weightSum;
}

}