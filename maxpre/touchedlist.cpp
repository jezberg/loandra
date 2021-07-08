#include <vector>
#include <deque>
#include <cstdint>
#include <map>
#include <iostream>

#include "touchedlist.hpp"
#include "global.hpp"
#include "probleminstance.hpp"

#define F first
#define S second

using namespace std;
namespace maxPreprocessor {

TouchedList::TouchedList(const ProblemInstance& pi_) : pi(pi_) {}

void TouchedList::init(int vars_) {
	vars = vars_;
	itr = 1;
	frontItr = 1;
	getItr = 1;
	touchedLiterals.clear();
	modLiterals.clear();
	touchedClauses.clear();
	modClauses.clear();
	techniqueItr.clear();
	getI.resize(vars*2);
	for (int i = 0; i < (int)getI.size(); i++) {
		getI[i] = 0;
	}
}

void TouchedList::shrink(bool log) {
	uint64_t mi = 1e18;
	for (auto t : techniqueItr) {
		mi = min(mi, t.S);
	}
	if (mi <= frontItr) return;
	if (mi - frontItr < (uint64_t)vars) return;
	frontItr = mi;
	int rm = 0;
	while (!touchedLiterals.empty() && touchedLiterals.front().F < frontItr) {
		touchedLiterals.pop_front();
		rm++;
	}
	touchedLiterals.shrink_to_fit();
	while (!touchedClauses.empty() && touchedClauses.front().F < frontItr) {
		touchedClauses.pop_front();
		rm++;
	}
	touchedClauses.shrink_to_fit();
	while (!modLiterals.empty() && modLiterals.front().F < frontItr) {
		modLiterals.pop_front();
		rm++;
	}
	modLiterals.shrink_to_fit();
	while (!modClauses.empty() && modClauses.front().F < frontItr) {
		modClauses.pop_front();
		rm++;
	}
	modClauses.shrink_to_fit();
	if (log) cerr<<"shrinked "<<rm<<endl;
}

void TouchedList::addVar() {
	vars++;
	getI.resize(vars*2);
}

void TouchedList::modClause(int c) {
	modClauses.push_back({itr++, c});
}

void TouchedList::modLiteral(int l) {
	modLiterals.push_back({itr++, l});
}

void TouchedList::touchClause(int c) {
	touchedClauses.push_back({itr++, c});
}

void TouchedList::touchLiteral(int l) {
	touchedLiterals.push_back({itr++, l});
}

void TouchedList::getTouchedLiteralsCh(uint64_t eitr, std::vector<int>& ret) {
	for (int i = (int)touchedLiterals.size() - 1; i >= 0 && touchedLiterals[i].F >= eitr; i--) {
		if (getI[touchedLiterals[i].S] != getItr) {
			getI[touchedLiterals[i].S] = getItr;
			ret.push_back(touchedLiterals[i].S);
		}
	}
	for (int i = (int)touchedClauses.size() - 1; i >= 0 && touchedClauses[i].F >= eitr; i--) {
		for (int l : pi.clauses[touchedClauses[i].S].lit) {
			if (getI[l] != getItr) {
				getI[l] = getItr;
				ret.push_back(l);
			}
		}
	}
}

void TouchedList::getModLiteralsCh(uint64_t eitr, std::vector<int>& ret) {
	for (int i = (int)modLiterals.size() - 1; i >= 0 && modLiterals[i].F >= eitr; i--) {
		if (getI[modLiterals[i].S] != getItr) {
			getI[modLiterals[i].S] = getItr;
			ret.push_back(modLiterals[i].S);
		}
	}
	for (int i = (int)modClauses.size() - 1; i >= 0 && modClauses[i].F >= eitr; i--) {
		for (int l : pi.clauses[modClauses[i].S].lit) {
			if (getI[l] != getItr) {
				getI[l] = getItr;
				ret.push_back(l);
			}
		}
	}
}

void TouchedList::getTouchedVarsCh(uint64_t eitr, std::vector<int>& ret) {
	for (int i = (int)touchedLiterals.size() - 1; i >= 0 && touchedLiterals[i].F >= eitr; i--) {
		if (getI[litVariable(touchedLiterals[i].S)] != getItr) {
			getI[litVariable(touchedLiterals[i].S)] = getItr;
			ret.push_back(litVariable(touchedLiterals[i].S));
		}
	}
	for (int i = (int)touchedClauses.size() - 1; i >= 0 && touchedClauses[i].F >= eitr; i--) {
		for (int l : pi.clauses[touchedClauses[i].S].lit) {
			if (getI[litVariable(l)] != getItr) {
				getI[litVariable(l)] = getItr;
				ret.push_back(litVariable(l));
			}
		}
	}
}

void TouchedList::getModVarsCh(uint64_t eitr, std::vector<int>& ret) {
	for (int i = (int)modLiterals.size() - 1; i >= 0 && modLiterals[i].F >= eitr; i--) {
		if (getI[litVariable(modLiterals[i].S)] != getItr) {
			getI[litVariable(modLiterals[i].S)] = getItr;
			ret.push_back(litVariable(modLiterals[i].S));
		}
	}
	for (int i = (int)modClauses.size() - 1; i >= 0 && modClauses[i].F >= eitr; i--) {
		for (int l : pi.clauses[modClauses[i].S].lit) {
			if (getI[litVariable(l)] != getItr) {
				getI[litVariable(l)] = getItr;
				ret.push_back(litVariable(l));
			}
		}
	}
}

vector<int> TouchedList::getModLiterals(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret(2*vars);
		for (int i = 0; i < 2*vars; i++) {
			ret[i] = i;
		}
		return ret;
	}
	else {
		getItr++;
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		getModLiteralsCh(eitr, ret);
		return ret;
	}
}

vector<int> TouchedList::getModVariables(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret(vars);
		for (int i = 0; i < vars; i++) {
			ret[i] = i;
		}
		return ret;
	}
	else {
		getItr++;
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		getModVarsCh(eitr, ret);
		return ret;
	}
}

vector<int> TouchedList::getTouchedLiterals(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret(2*vars);
		for (int i = 0; i < 2*vars; i++) {
			ret[i] = i;
		}
		return ret;
	}
	else {
		getItr++;
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		getModLiteralsCh(eitr, ret);
		getTouchedLiteralsCh(eitr, ret);
		return ret;
	}
}

vector<int> TouchedList::getTouchedVariables(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret(vars);
		for (int i = 0; i < vars; i++) {
			ret[i] = i;
		}
		return ret;
	}
	else {
		getItr++;
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		getModVarsCh(eitr, ret);
		getTouchedVarsCh(eitr, ret);
		return ret;
	}
}

vector<int> TouchedList::getBinaryLiterals(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret(2*vars);
		for (int i = 0; i < 2*vars; i++) {
			ret[i] = i;
		}
		return ret;
	}
	else {
		getItr++;
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		for (int i = (int)modClauses.size() - 1; i >= 0 && modClauses[i].F >= eitr; i--) {
			if (pi.clauses[modClauses[i].S].lit.size() == 2) {
				for (int l : pi.clauses[modClauses[i].S].lit) {
					if (getI[l] != getItr) {
						getI[l] = getItr;
						ret.push_back(l);
					}
				}
			}
		}
		return ret;
	}
}

vector<int> TouchedList::getModClauses(string technique) {
	if (techniqueItr.count(technique) == 0) {
		techniqueItr[technique] = itr;
		vector<int> ret;
		ret.reserve(pi.clauses.size());
		for (int i = 0; i < (int)pi.clauses.size(); i++) {
			if (!pi.isClauseRemoved(i)) {
				ret.push_back(i);
			}
		}
		return ret;
	}
	else {
		getItr++;
		if (getI.size() < pi.clauses.size()) {
			getI.resize(pi.clauses.size());
		}
		uint64_t eitr = techniqueItr[technique];
		techniqueItr[technique] = itr;
		vector<int> ret;
		for (int i = (int)modClauses.size() - 1; i >= 0 && modClauses[i].F >= eitr; i--) {
			if (!pi.isClauseRemoved(modClauses[i].S)) {
				if (getI[modClauses[i].S] != getItr) {
					getI[modClauses[i].S] = getItr;
					ret.push_back(modClauses[i].S);
				}
			}
		}
		return ret;
	}
}

void TouchedList::setItr(string technique) {
	techniqueItr[technique] = itr;
}

}