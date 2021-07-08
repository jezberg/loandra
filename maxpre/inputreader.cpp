#include <iostream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <limits>
#include <cassert>

#include "inputreader.hpp"
#include "global.hpp"
#include "cardinalityconstraint.hpp"

using namespace std;
namespace maxPreprocessor {

int InputReader::readWeightedClause(istream& input) {
	uint64_t weight;
	vector<int> clause;
	input>>weight;
	if (!input.good()) {
		return 1;
	}
	if (weight > top) weight = top;
	int lit;
	while (input>>lit) {
		if (lit == 0) break;
		clause.push_back(lit);
	}
	if (input.fail()) {
		return 1;
	}
	if (weight == 0) return 0;
	clauses.push_back(clause);
	weights.push_back(weight);
	return 0;
}

int InputReader::readClause(istream& input, uint64_t defaultWeight) {
	vector<int> clause;
	int lit;
	while (input>>lit) {
		if (lit == 0) break;
		clause.push_back(lit);
	}
	if (input.fail()) {
		return 1;
	}
	clauses.push_back(clause);
	weights.push_back(defaultWeight);
	return 0;
}

int InputReader::readCardinalityConstraint(istream& inputF) {
	string inputLine;
	getline(inputF, inputLine);
	stringstream input(inputLine);
	string coding;
	input>>coding;
	if (input.fail()) return 1;
	int codingN;
	if (coding == "CARD") {
		codingN = CardinalityConstraint::CODING_DEFAULT;
	}
	else {
		return 1;
	}
	
	vector<int> inputLiterals;
	vector<int> outputLiterals;
	vector<int> addLiterals;
	string in;
	string compr;
	int p = 0;
	int k = 0;
	bool foundK = false;
	bool iff = 0;
	while (!input.eof() && input>>in) {
		if (input.fail()) return 1;
		if (in == "=" || in == ">=" || in == "<=" || in == ">" || in == "<" || in == "==" || in == "!=" || in == "<:" || in == ">:" || in == "::") {
			if (p != 0) return 1;
			compr = in;
			p = 1;
		}
		else if (in == "OUT") {
			if (p > 1) return 1;
			p = 2;
		}
		else if (in == "ADD") {
			if (p == 0 || p > 2) return 1;
			p = 3;
		}
		else if (in == "IFF") {
			iff = true;
		}
		else {
			stringstream ss;
			ss<<in;
			int t;
			ss>>t;
			if (!ss.eof()) return 1;
			if (ss.fail()) return 1;
			if (p == 0) {
				if (t == 0) return 1;
				inputLiterals.push_back(t);
			}
			else if (p == 1) {
				if (foundK && t == 0) break;
				if (foundK) return 1;
				foundK = true;
				k = t;
			}
			else if (p == 2) {
				if (t == 0) break;
				outputLiterals.push_back(t);
			}
			else if (p == 3) {
				if (t == 0) break;
				addLiterals.push_back(t);
			}
		}
		input>>ws;
	}
	if (!foundK) return 1;
	if (input.fail()) return 1;
	if (inputLiterals.size() == 0) return 1;
	if (compr == ">") {
		compr = ">=";
		k++;
	}
	if (compr == "<") {
		compr = "<=";
		k--;
	}
	if (compr == "==") compr = "=";
	
	int type = -1;
	
	if (compr == ">=") {
		compr = "<=";
		k = (int)inputLiterals.size() - k;
		for (int& lit : inputLiterals) {
			lit = -lit;
		}
	}
	
	if (compr == "<=") type = CardinalityConstraint::CONSTRAINT_LESS_EQ;
	else if (compr == "=") type = CardinalityConstraint::CONSTRAINT_EQUAL;
	else if (compr == "!=") type = CardinalityConstraint::CONSTRAINT_NEQUAL;
	else if (compr == "<:") type = CardinalityConstraint::CONSTRAINT_ALL_LESS;
	else if (compr == ">:") type = CardinalityConstraint::CONSTRAINT_ALL_GREATER;
	else if (compr == "::") type = CardinalityConstraint::CONSTRAINT_ALL_BOTH;
	else return 1;
	
	if (type == CardinalityConstraint::CONSTRAINT_ALL_LESS || type == CardinalityConstraint::CONSTRAINT_ALL_GREATER || type == CardinalityConstraint::CONSTRAINT_ALL_BOTH) {
		if (k != (int)outputLiterals.size()) return 1;
		if (k > (int)inputLiterals.size()) return 1;
		if (k <= 0) return 1;
	}
	else {
		if (outputLiterals.size() > 1) return 1;
	}
	
	if (outputLiterals.size() == 0) {
		iff = false;
	}
	if (type == CardinalityConstraint::CONSTRAINT_ALL_BOTH) {
		iff = true;
	}
	cardinalityConstraints.push_back(CardinalityConstraint(type, codingN, inputLiterals, outputLiterals, k, addLiterals, iff));
	if (input.fail()) return 1;
	return 0;
}

int InputReader::readLine(istream& input) {
	input>>ws;
	while (input.peek() == 'c') {
		currentLine++;
		input.ignore(numeric_limits<streamsize>::max(), '\n');
		input>>ws;
	}
	if (input.eof()) return 1;
	int status;
	if (input.peek() == 'C') {
		status = readCardinalityConstraint(input);
		if (status == 1) readError = "Failed to read a cardinality constraint (line "+to_string(currentLine)+")";
	}
	else if (inputFormat == INPUT_FORMAT_WPMS) {
		status = readWeightedClause(input);
		if (status == 1) readError = "Failed to read a clause (line "+to_string(currentLine)+")";
	}
	else if (inputFormat == INPUT_FORMAT_MS) {
		status = readClause(input, 1);
		if (status == 1) readError = "Failed to read a clause (line "+to_string(currentLine)+")";
	}
	else if (inputFormat == INPUT_FORMAT_SAT) {
		status = readClause(input, top);
		if (status == 1) readError = "Failed to read a clause (line "+to_string(currentLine)+")";
	}
	else {
		status = 1;
		readError = "Wrong input format (line "+to_string(currentLine)+")";
	}
	assert(status == 0 || status == 1);
	if (status == 0) return 0;
	else return 2;
}

int InputReader::readClauses(istream& input, bool maxSat) {
	clauses.clear();
	weights.clear();
	cardinalityConstraints.clear();
	readError = "";
	currentLine = 0;
	string line;
	while (getline(input, line)) {
		currentLine++;
		if (line.size() > 0 && line[0] == 'p') break;
	}
	if (line.size() == 0 || line[0] != 'p') {
		readError = "Cannot find the p-line of the input";
		return 1;
	}
	stringstream parameterLine;
	parameterLine<<line;
	parameterLine>>line;
	if (line != "p") {
		readError = "Invalid p-line of the input (line "+to_string(currentLine)+")";
		return 1;
	}
	parameterLine>>line>>nbvar>>nbclauses;
	clauses.reserve(nbclauses);
	weights.reserve(nbclauses);
	if (parameterLine.fail()) {
		readError = "Invalid p-line of the input (line "+to_string(currentLine)+")";
		return 1;
	}
	if (line == "wcnf" && maxSat) {
		inputFormat = INPUT_FORMAT_WPMS;
		if (!(parameterLine>>top)) {
			cerr<<"Weighted Max-SAT input format"<<endl;
			top = (1ull << 63ull) - 1;
		}
		else {
			cerr<<"Weighted Partial Max-SAT input format"<<endl;
		}
		cerr<<"Top "<<top<<endl;
	}
	else if (line == "cnf" && maxSat) {
		inputFormat = INPUT_FORMAT_MS;
		cerr<<"Max-SAT input format"<<endl;
		top = (1ull << 63ull) - 1;
	}
	else if (line == "cnf" && !maxSat) {
		inputFormat = INPUT_FORMAT_SAT;
		cerr<<"SAT input format"<<endl;
		top = (1ull << 63ull) - 1;
	}
	else {
		readError = "Invalid p-line of the input (problem type "+line+" and the dimacs format do not match) (line "+to_string(currentLine)+")";
		return 1;
	}
	while (1) {
		currentLine++;
		int status = readLine(input);
		if (status == 1) break;
		else if (status == 2) return 1;
		assert(status == 0);
	}
	int freeVar = 1;
	for (auto& clause : clauses) {
		for (int lit : clause) {
			freeVar = max(freeVar, abs(lit) + 1);
		}
	}
	for (auto& cc : cardinalityConstraints) {
		for (int lit : cc.inputLiterals) freeVar = max(freeVar, abs(lit) + 1);
		for (int lit : cc.outputLiterals) freeVar = max(freeVar, abs(lit) + 1);
		for (int lit : cc.addLiterals) freeVar = max(freeVar, abs(lit) + 1);
	}
	for (auto& cc : cardinalityConstraints) {
		vector<vector<int> > cls = cc.getClauses(freeVar);
		for (auto& clause : cls) {
			clauses.push_back(clause);
			weights.push_back(top);
			for (int lit : clause) {
				freeVar = max(freeVar, abs(lit) + 1);
			}
		}
	}
	return 0;
}

}