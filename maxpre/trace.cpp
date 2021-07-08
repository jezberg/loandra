#include <vector>
#include <iostream>
#include <cstdint>
#include <cassert>

#include "trace.hpp"
#include "utility.hpp"
#include "global.hpp"

#define F first
#define S second

using namespace std;
namespace maxPreprocessor {

void Trace::setVar(int var, bool value) {
	operations.push_back(1);
	data.push_back({var, value});
}

// The variable and the clauses which contained negation of var
void Trace::BVE(int var, const vector<vector<int> >& nClauses) {
	operations.push_back(2);
	data.push_back(vector<int>());
	data.back().push_back(var);
	for (auto& c : nClauses) {
		for (int lit : c) {
			if (litVariable(lit) != var) data.back().push_back(lit);
		}
		data.back().push_back(-1);
	}
}

void Trace::BCE(int lit, const vector<int>& clause) {
	operations.push_back(3);
	data.push_back(vector<int>());
	data.back().push_back(lit);
	for (int l : clause) {
		if (l != lit) {
			data.back().push_back(l);
		}
	}
}

void Trace::BCR(int lbl1, int lbl2, const vector<vector<int> >& nClauses) {
	operations.push_back(4);
	data.push_back(vector<int>());
	data.back().push_back(lbl1);
	data.back().push_back(lbl2);
	for (auto& c : nClauses) {
		for (int lit : c) {
			assert(litVariable(lit) != litVariable(lbl2));
			if (litVariable(lit) != litVariable(lbl1)) data.back().push_back(lit);
		}
		data.back().push_back(-1);
	}
}
// lit2 = lit1
void Trace::setEqual(int lit1, int lit2) {
	operations.push_back(5);
	data.push_back({lit1, lit2});
}

void Trace::LS(int lbl, int lit, const vector<int>& clause) {
	operations.push_back(6);
	data.push_back(vector<int>());
	data.back().push_back(lbl);
	data.back().push_back(lit);
	for (int l : clause) {
		if (l != lit) {
			data.back().push_back(l);
		}
	}
}

void Trace::labelEliminate(int lbl1, int lbl2, int tautli) {
	operations.push_back(7);
	data.push_back({lbl1, lbl2, tautli});
}

void Trace::removeWeight(uint64_t weight) {
	if (weight != HARDWEIGHT) {
		removedWeight += weight;
	}
}

pair<vector<int>, uint64_t> Trace::getSolution(const vector<int>& trueLits, uint64_t weight, int vars, int originalVars) {
	vector<int> value(vars);
	for (int lit : trueLits) {
		if (abs(lit) <= vars) {
			if (lit < 0) value[(-lit) - 1] = VAR_FALSE;
			else value[lit - 1] = VAR_TRUE;
		}
	}
	
	for (int i = 0; i < vars; i++) {
		if (value[i] == VAR_UNDEFINED) {
			value[i] = VAR_TRUE;
		}
	}
	
	for (int i = (int)operations.size() -1; i >= 0; i--) {
		if (operations[i] == 1) {
			assert(data[i][0] < vars);
			
			if (data[i][1] == true) {
				value[data[i][0]] = VAR_TRUE;
			}
			else {
				value[data[i][0]] = VAR_FALSE;
			}
		}
		else if(operations[i] == 2) {
			int var = data[i][0];
			bool f = false;
			bool ff = false;
			for (int j = 1; j < (int)data[i].size(); j++) {
				if (data[i][j] == -1) {
					if (!f) {
						value[var] = VAR_FALSE;
						ff = true;
						break;
					}
					f = false;
				}
				else {
					if (value[litVariable(data[i][j])] == VAR_UNDEFINED) {
						if (litValue(data[i][j]) == true) value[litVariable(data[i][j])] = VAR_TRUE;
						else value[litVariable(data[i][j])] = VAR_FALSE;
					}
					if ((litValue(data[i][j]) == true && value[litVariable(data[i][j])] == VAR_TRUE) || 
						(litValue(data[i][j]) == false && value[litVariable(data[i][j])] == VAR_FALSE)) {
						f = true;
					}
				}
			}
			if (!ff) value[var] = VAR_TRUE;
		}
		else if(operations[i] == 3) {
			bool sat = false;
			for (int j = 1; j < (int)data[i].size(); j++) {
				if (value[litVariable(data[i][j])] == VAR_UNDEFINED) {
					if (litValue(data[i][j]) == true) value[litVariable(data[i][j])] = VAR_TRUE;
					else value[litVariable(data[i][j])] = VAR_FALSE;
				}
				if ((litValue(data[i][j]) == true && value[litVariable(data[i][j])] == VAR_TRUE) || 
					(litValue(data[i][j]) == false && value[litVariable(data[i][j])] == VAR_FALSE)) {
					sat = true;
					break;
				}
			}
			if (!sat) {
				if (litValue(data[i][0]) == false) value[litVariable(data[i][0])] = VAR_FALSE;
				else value[litVariable(data[i][0])] = VAR_TRUE;
			}
		}
		else if(operations[i] == 4) {
			bool f = false;
			bool ff = false;
			int lbl1 = data[i][0];
			int lbl2 = data[i][1];
			for (int j = 2; j < (int)data[i].size(); j++) {
				if (data[i][j] == -1) {
					if (!f) {
						if (litValue(lbl1) == true) value[litVariable(lbl1)] = VAR_TRUE;
						else value[litVariable(lbl1)] = VAR_FALSE;
						ff = true;
						break;
					}
					f = false;
				}
				else {
					if (value[litVariable(data[i][j])] == VAR_UNDEFINED) {
						if (litValue(data[i][j]) == true) value[litVariable(data[i][j])] = VAR_TRUE;
						else value[litVariable(data[i][j])] = VAR_FALSE;
					}
					if ((litValue(data[i][j]) == true && value[litVariable(data[i][j])] == VAR_TRUE) || 
						(litValue(data[i][j]) == false && value[litVariable(data[i][j])] == VAR_FALSE)) {
						f = true;
					}
				}
			}
			if (!ff) {
				if (litValue(lbl1) == true) value[litVariable(lbl1)] = VAR_FALSE;
				else value[litVariable(lbl1)] = VAR_TRUE;
				if (litValue(lbl2)== true) value[litVariable(lbl2)] = VAR_TRUE;
				else value[litVariable(lbl2)] = VAR_FALSE;
			}
		}
		else if(operations[i] == 5) {
			if (value[litVariable(data[i][0])] == VAR_UNDEFINED) {
				value[litVariable(data[i][0])] = VAR_TRUE;
			}
			if (litValue(data[i][0]) != litValue(data[i][1])) {
				if (value[litVariable(data[i][0])] == VAR_TRUE) {
					value[litVariable(data[i][1])] = VAR_FALSE;
				}
				else {
					value[litVariable(data[i][1])] = VAR_TRUE;
				}
			}
			else {
				if (value[litVariable(data[i][0])] == VAR_TRUE) {
					value[litVariable(data[i][1])] = VAR_TRUE;
				}
				else {
					value[litVariable(data[i][1])] = VAR_FALSE;
				}
			}
		}
		else if(operations[i] == 6) {
			if ((litValue(data[i][0]) == true && value[litVariable(data[i][0])] == VAR_TRUE) ||
				(litValue(data[i][0]) == false && value[litVariable(data[i][0])] == VAR_FALSE)) {
				bool sat = false;
				for (int j = 2; j < (int)data[i].size(); j++) {
					if (value[litVariable(data[i][j])] == VAR_UNDEFINED) {
						if (litValue(data[i][j]) == true) value[litVariable(data[i][j])] = VAR_TRUE;
						else value[litVariable(data[i][j])] = VAR_FALSE;
					}
					if ((litValue(data[i][j]) == true && value[litVariable(data[i][j])] == VAR_TRUE) || 
						(litValue(data[i][j]) == false && value[litVariable(data[i][j])] == VAR_FALSE)) {
						sat = true;
						break;
					}
				}
				if (!sat) {
					if (litValue(data[i][1]) == false) value[litVariable(data[i][1])] = VAR_FALSE;
					else value[litVariable(data[i][1])] = VAR_TRUE;
				}
			}
		}
		else if(operations[i] == 7) {
			if ((litValue(data[i][0]) == true && value[litVariable(data[i][0])] == VAR_FALSE) ||
				(litValue(data[i][0]) == false && value[litVariable(data[i][0])] == VAR_TRUE)) {
				if ((litValue(data[i][2]) == true && value[litVariable(data[i][2])] == VAR_TRUE) ||
				(litValue(data[i][2]) == false && value[litVariable(data[i][2])] == VAR_FALSE)) {
					if (litValue(data[i][1]) == false) {
						value[litVariable(data[i][1])] = VAR_FALSE;
					}
					else {
						value[litVariable(data[i][1])] = VAR_TRUE;
					}
					if (litValue(data[i][0]) == false) {
						value[litVariable(data[i][0])] = VAR_TRUE;
					}
					else {
						value[litVariable(data[i][0])] = VAR_FALSE;
					}
				}
				else {
					if (litValue(data[i][1]) == false) {
						value[litVariable(data[i][1])] = VAR_TRUE;
					}
					else {
						value[litVariable(data[i][1])] = VAR_FALSE;
					}
					if (litValue(data[i][0]) == false) {
						value[litVariable(data[i][0])] = VAR_FALSE;
					}
					else {
						value[litVariable(data[i][0])] = VAR_TRUE;
					}
				}
			}
		}
		else {
			assert(0);
		}
	}
	
	vector<int> retLit;
	for (int i = 0; i < originalVars; i++) {
		if (value[i] == VAR_FALSE) {
			retLit.push_back(-(i + 1));
		}
		else if(value[i] == VAR_TRUE) {
			retLit.push_back(i + 1);
		}
		else {
			retLit.push_back(i + 1);
		}
	}
	return {retLit, weight};
}

void Trace::printSolution(ostream& output, const vector<int>& trueLits, uint64_t weight, int vars, int originalVars) {
	auto solution = getSolution(trueLits, weight, vars, originalVars);
	output << "v ";
	for (int lit : solution.F) {
		output << lit << " ";
	}
	output << '\n';
	output << "s OPTIMUM FOUND\n";
	output << "o " << solution.S << '\n';
	
	output.flush();
}

}