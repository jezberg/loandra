#include <vector>
#include <cassert>
#include <iostream>
#include <map>

#define F first
#define S second

#include "cardinalityconstraint.hpp"

using namespace std;
namespace maxPreprocessor {
CardinalityConstraint::CardinalityConstraint(int type_, int coding_, const vector<int>& inputLiterals_, const vector<int>& outputLiterals_, int K_, const vector<int>& addLiterals_, bool iff_) :
type(type_), coding(coding_), inputLiterals(inputLiterals_), outputLiterals(outputLiterals_), K(K_), addLiterals(addLiterals_), iff(iff_ || type == CONSTRAINT_ALL_BOTH) {
	assert(type == CONSTRAINT_LESS_EQ || type == CONSTRAINT_EQUAL || type == CONSTRAINT_ALL_LESS || type == CONSTRAINT_NEQUAL || type == CONSTRAINT_ALL_GREATER || type == CONSTRAINT_ALL_BOTH);
	assert(coding == CODING_NETWORK);
	assert((int)inputLiterals.size() > 0);
}
CardinalityConstraint::Formula::Formula(int freeVar_) : freeVar(freeVar_) {}
size_t CardinalityConstraint::Formula::size() const {
	return clauses.size();
}
vector<int> CardinalityConstraint::HMerge(const vector<int>& A, const vector<int>& B, Formula& F, bool d1, bool d2) {
	// use the bit reversal trick to do HMerge without recursion
	assert(A.size() == B.size());
	int n = A.size();
	vector<int> c(2*n);
	for (int i = 0; i < n; i++) {
		c[2*i] = F.freeVar++;
		c[2*i+1] = F.freeVar++;
		if (d1) {
			F.clauses.push_back({-A[i], -B[i], c[2*i+1]});
			F.clauses.push_back({-A[i], c[2*i]});
			F.clauses.push_back({-B[i], c[2*i]});
		}
		if (d2) {
			F.clauses.push_back({A[i], B[i], -c[2*i]});
			F.clauses.push_back({A[i], -c[2*i+1]});
			F.clauses.push_back({B[i], -c[2*i+1]});
		}
	}
	for (int i = 0; i < n; i++) {
		int u = 0;
		for (int j = 1; j < n; j *= 2) {
			u *= 2;
			if (i & j) u++;
		}
		if (i < u) {
			swap(c[2*i], c[2*u]);
			swap(c[2*i+1], c[2*u+1]);
		}
	}
	for (int le = 4; le <= 2*n; le *= 2) {
		for (int i = 0; i < 2*n; i+=le) {
			for (int j = 0; j+1 < le/2; j++) {
				if (d1) {
					F.clauses.push_back({-c[i+j+1], -c[i+j+le/2], F.freeVar+2*j+1});
					F.clauses.push_back({-c[i+j+1], F.freeVar+2*j});
					F.clauses.push_back({-c[i+j+le/2], F.freeVar+2*j});
				}
				if (d2) {
					F.clauses.push_back({c[i+j+1], c[i+j+le/2], -(F.freeVar+2*j)});
					F.clauses.push_back({c[i+j+1], -(F.freeVar+2*j+1)});
					F.clauses.push_back({c[i+j+le/2], -(F.freeVar+2*j+1)});
				}
			}
			for (int j = 1; j+1 < le; j++) c[i+j] = F.freeVar++;
		}
	}
	return c;
}
vector<int> CardinalityConstraint::SMerge(const vector<int>& A, const vector<int>& B, Formula& F, bool d1, bool d2) {
	// use the bit reversal trick to do SMerge without recursion
	assert(A.size() == B.size());
	int n = A.size();
	vector<int> c(2*n);
	for (int i = 0; i < n; i++) {
		c[2*i] = F.freeVar++;
		if (i < n/2) c[2*i+1] = F.freeVar++;
		if (d1) {
			if (i < n/2) F.clauses.push_back({-A[i], -B[i], c[2*i+1]});
			F.clauses.push_back({-A[i], c[2*i]});
			F.clauses.push_back({-B[i], c[2*i]});
		}
		if (d2) {
			F.clauses.push_back({A[i], B[i], -c[2*i]});
			if (i < n/2) F.clauses.push_back({A[i], -c[2*i+1]});
			if (i < n/2) F.clauses.push_back({B[i], -c[2*i+1]});
		}
	}
	for (int i = 0; i < n; i++) {
		int u = 0;
		for (int j = 1; j < n; j *= 2) {
			u *= 2;
			if (i & j) u++;
		}
		if (i < u) {
			swap(c[2*i], c[2*u]);
			swap(c[2*i+1], c[2*u+1]);
		}
	}
	int oc = 0;
	for (int le = 4; le <= 2*n; le *= 2) {
		for (int i = 0; i < 2*n; i+=le) {
			if (le == 2*n) oc = 1;
			for (int j = 0; j < le/4; j++) {
				if (d1) {
					if (j+1 < le/4 || !oc) F.clauses.push_back({-c[i+j+1], -c[i+j+le/2], F.freeVar+2*j+1});
					F.clauses.push_back({-c[i+j+1], F.freeVar+2*j});
					F.clauses.push_back({-c[i+j+le/2], F.freeVar+2*j});
				}
				if (d2) {
					F.clauses.push_back({c[i+j+1], c[i+j+le/2], -(F.freeVar+2*j)});
					if (j+1 < le/4 || !oc) F.clauses.push_back({c[i+j+1], -(F.freeVar+2*j+1)});
					if (j+1 < le/4 || !oc) F.clauses.push_back({c[i+j+le/2], -(F.freeVar+2*j+1)});
				}
			}
			for (int j = 1; j < le/2; j++) c[i+j] = F.freeVar++;
			if (!oc) c[i+le/2] = F.freeVar++;
			oc ^= 1;
		}
	}
	c.resize(n);
	return c;
}
vector<int> CardinalityConstraint::HSort(const vector<int>& A, int l, int r, Formula& F, bool d1, bool d2) {
	int le = r-l+1;
	assert(le >= 2);
	assert(le%2 == 0);
	if (le == 2) {
		vector<int> a = {A[l]};
		vector<int> b = {A[r]};
		return HMerge(a, b, F, d1, d2);
	}
	else {
		int m = (l+r)/2;
		vector<int> a = HSort(A, l, m, F, d1, d2);
		vector<int> b = HSort(A, m+1, r, F, d1, d2);
		return HMerge(a, b, F, d1, d2);
	}
}
vector<int> CardinalityConstraint::Card(const vector<int>& A, int l, int r, int k, Formula& F, bool d1, bool d2) {
	int n = (r-l+1);
	assert(n%k == 0);
	if (n == k) {
		return HSort(A, l, r, F, d1, d2);
	}
	else {
		vector<int> a = Card(A, l, l+k-1, k, F, d1, d2);
		vector<int> b = Card(A, l+k, r, k, F, d1, d2);
		assert((int)a.size() == k && (int)b.size() == k);
		return SMerge(a, b, F, d1, d2);
	}
}

pair<CardinalityConstraint::Formula, int> CardinalityConstraint::getLessNetwork(int freeVar, vector<int> input) {
	int k = 1;
	while (k <= K) k *= 2;
	Formula F(freeVar);
	while ((int)input.size()%k != 0) {
		F.clauses.push_back({-F.freeVar});
		input.push_back(F.freeVar++);
	}
	vector<int> output = Card(input, 0, (int)input.size()-1, k, F, true, iff);
	return {F, -output[K]};
}

pair<CardinalityConstraint::Formula, int> CardinalityConstraint::getGreaterNetwork(int freeVar, vector<int> input) {
	int KK = (int)input.size() - K;
	int k = 1;
	while (k < KK) k *= 2;
	Formula F(freeVar);
	for (int& lit : input) {
		lit = -lit;
	}
	while ((int)input.size()%k != 0) {
		F.clauses.push_back({-F.freeVar});
		input.push_back(F.freeVar++);
	}
	vector<int> output = Card(input, 0, (int)input.size()-1, k, F, iff, true);
	return {F, output[KK-1]};
}
pair<CardinalityConstraint::Formula, int> CardinalityConstraint::getEqualNetwork(int freeVar, vector<int> input) {
	int KK = K;
	if ((int)input.size() - K < KK) {
		KK = (int)input.size() - K;
		for (int& lit : input) {
			lit = -lit;
		}
	}
	int k = 1;
	while (k <= KK) k *= 2;
	Formula F(freeVar);
	while ((int)input.size()%k != 0) {
		F.clauses.push_back({-F.freeVar});
		input.push_back(F.freeVar++);
	}
	vector<int> output = Card(input, 0, (int)input.size()-1, k, F, true, true);
	int oV = F.freeVar++;
	F.clauses.push_back({output[KK-1], -oV});
	F.clauses.push_back({-output[KK], -oV});
	if (iff) {
		F.clauses.push_back({-output[KK-1], output[KK], oV});
	}
	return {F, oV};
}
pair<CardinalityConstraint::Formula, int> CardinalityConstraint::getNequalNetwork(int freeVar, vector<int> input) {
	int KK = K;
	if ((int)input.size() - K < KK) {
		KK = (int)input.size() - K;
		for (int& lit : input) {
			lit = -lit;
		}
	}
	int k = 1;
	while (k <= KK) k *= 2;
	Formula F(freeVar);
	while ((int)input.size()%k != 0) {
		F.clauses.push_back({-F.freeVar});
		input.push_back(F.freeVar++);
	}
	vector<int> output = Card(input, 0, (int)input.size()-1, k, F, true, true);
	int oV = F.freeVar++;
	F.clauses.push_back({-output[KK-1], output[KK], oV});
	if (iff) {
		F.clauses.push_back({output[KK-1], -oV});
		F.clauses.push_back({-output[KK], -oV});
	}
	return {F, -oV};
}
pair<CardinalityConstraint::Formula, int> CardinalityConstraint::getBestLessNetwork(int freeVar, vector<int> input) {
	auto cLess = getLessNetwork(freeVar, input);
	auto cGreater = getGreaterNetwork(freeVar, input);
	if (cLess.F.size() <= cGreater.F.size()) {
		return cLess;
	}
	else {
		return cGreater;
	}
}
vector<vector<int> > CardinalityConstraint::getClauses(int freeVar) {
	// handle all special/easy cases separately
	int n = inputLiterals.size();
	assert(n > 0);
	if (type == CONSTRAINT_NEQUAL) {
		if (K < 0 || K > n) {
			if (outputLiterals.size() == 0) {
				return {};
			}
			else if (outputLiterals.size() == 1) {
				if (iff) return {{outputLiterals[0]}};
				else return {};
			}
			else {
				assert(0);
			}
		}
		if (K == 0) {
			if (outputLiterals.size() == 0) {
				return {inputLiterals};
			}
			else if (outputLiterals.size() == 1) {
				auto retC = inputLiterals;
				retC.push_back(-outputLiterals[0]);
				if (iff) {
					vector<vector<int> > clauses = {retC};
					for (int lit : inputLiterals) {
						clauses.push_back({-lit, outputLiterals[0]});
					}
					return clauses;
				}
				else {
					return {retC};
				}
			}
			else {
				assert(0);
			}
		}
		if (K == n) {
			if (outputLiterals.size() == 0) {
				auto retC = inputLiterals;
				for (int& lit : retC) {
					lit = -lit;
				}
				return {retC};
			}
			else if (outputLiterals.size() == 1) {
				auto retC = inputLiterals;
				for (int& lit : retC) {
					lit = -lit;
				}
				retC.push_back(-outputLiterals[0]);
				if (iff) {
					vector<vector<int> > clauses = {retC};
					for (int lit : inputLiterals) {
						clauses.push_back({lit, outputLiterals[0]});
					}
					return clauses;
				}
				else {
					return {retC};
				}
			}
			else {
				assert(0);
			}
		}
	}
	else if (type == CONSTRAINT_EQUAL) {
		if (K < 0 || K > n) {
			if (outputLiterals.size() == 0) {
				return {{freeVar}, {-freeVar}};
			}
			else if (outputLiterals.size() == 1) {
				return {{-outputLiterals[0]}};
			}
			else {
				assert(0);
			}
		}
		if (K == 0) {
			if (outputLiterals.size() == 0) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({-lit});
				}
				return clauses;
			}
			else if (outputLiterals.size() == 1) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({-lit, -outputLiterals[0]});
				}
				if (iff) {
					vector<int> retC = inputLiterals;
					retC.push_back(outputLiterals[0]);
					clauses.push_back(retC);
				}
				return clauses;
			}
			else {
				assert(0);
			}
		}
		if (K == n) {
			if (outputLiterals.size() == 0) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({lit});
				}
				return clauses;
			}
			else if (outputLiterals.size() == 1) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({lit, -outputLiterals[0]});
				}
				if (iff) {
					vector<int> retC = inputLiterals;
					for (int& lit : retC) {
						lit = -lit;
					}
					retC.push_back(outputLiterals[0]);
					clauses.push_back(retC);
				}
				return clauses;
			}
			else {
				assert(0);
			}
		}
	}
	else if (type == CONSTRAINT_LESS_EQ) {
		if (K < 0) {
			if (outputLiterals.size() == 0) {
				return {{freeVar}, {-freeVar}};
			}
			else if (outputLiterals.size() == 1) {
				return {{-outputLiterals[0]}};
			}
			else {
				assert(0);
			}
		}
		if (K >= n) {
			if (outputLiterals.size() == 0) {
				return {};
			}
			else if (outputLiterals.size() == 1) {
				if (iff) return {{outputLiterals[0]}};
				else return {};
			}
			else {
				assert(0);
			}
		}
		if (K == 0) {
			if (outputLiterals.size() == 0) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({-lit});
				}
				return clauses;
			}
			else if (outputLiterals.size() == 1) {
				vector<vector<int> > clauses;
				for (int lit : inputLiterals) {
					clauses.push_back({-lit, -outputLiterals[0]});
				}
				if (iff) {
					vector<int> retC = inputLiterals;
					retC.push_back(outputLiterals[0]);
					clauses.push_back(retC);
				}
				return clauses;
			}
			else {
				assert(0);
			}
		}
		if (K == n-1) {
			if (outputLiterals.size() == 0) {
				auto retC = inputLiterals;
				for (int& lit: retC) {
					lit = -lit;
				}
				return {retC};
			}
			else if (outputLiterals.size() == 1) {
				auto retC = inputLiterals;
				for (int& lit: retC) {
					lit = -lit;
				}
				retC.push_back(-outputLiterals[0]);
				if (iff) {
					vector<vector<int> > clauses = {retC};
					for (int lit : inputLiterals) {
						clauses.push_back({lit, outputLiterals[0]});
					}
					return clauses;
				}
				else {
					return {retC};
				}
			}
			else {
				assert(0);
			}
		}
	}
	else if (type == CONSTRAINT_ALL_LESS || type == CONSTRAINT_ALL_GREATER || type == CONSTRAINT_ALL_BOTH) {
		if (n == 1) {
			assert(outputLiterals.size() == 1);
			vector<vector<int> > clauses;
			if (type == CONSTRAINT_ALL_LESS || iff) {
				clauses.push_back({-inputLiterals[0], outputLiterals[0]});
			}
			if (type == CONSTRAINT_ALL_GREATER || iff) {
				clauses.push_back({inputLiterals[0], -outputLiterals[0]});
			}
			return clauses;
		}
	}
	else assert(0);
	assert((int)inputLiterals.size() > 1);
	// now n > 1
	vector<vector<int> > clauses;
	if (type == CONSTRAINT_LESS_EQ) {
		assert(0 < K && K+1 < (int)inputLiterals.size());
		if (coding == CODING_NETWORK) {
			pair<Formula, int> network = getBestLessNetwork(freeVar, inputLiterals);
			clauses = network.F.clauses;
			if (outputLiterals.size() == 0) {
				clauses.push_back({network.S});
			}
			else if (outputLiterals.size() == 1) {
				clauses.push_back({network.S, -outputLiterals[0]});
				if (iff) {
					clauses.push_back({-network.S, outputLiterals[0]});
				}
			}
			else {
				assert(0);
			}
		}
		else {
			assert(0);
		}
	}
	else if (type == CONSTRAINT_EQUAL) {
		assert(0 < K && K < (int)inputLiterals.size());
		if (coding == CODING_NETWORK) {
			pair<Formula, int> network = getEqualNetwork(freeVar, inputLiterals);
			clauses = network.F.clauses;
			if (outputLiterals.size() == 0) {
				clauses.push_back({network.S});
			}
			else if (outputLiterals.size() == 1) {
				clauses.push_back({network.S, -outputLiterals[0]});
				if (iff) {
					clauses.push_back({-network.S, outputLiterals[0]});
				}
			}
			else {
				assert(0);
			}
		}
		else {
			assert(0);
		}
	}
	else if (type == CONSTRAINT_NEQUAL) {
		assert(0 < K && K < (int)inputLiterals.size());
		if (coding == CODING_NETWORK) {
			pair<Formula, int> network = getNequalNetwork(freeVar, inputLiterals);
			clauses = network.F.clauses;
			if (outputLiterals.size() == 0) {
				clauses.push_back({network.S});
			}
			else if (outputLiterals.size() == 1) {
				clauses.push_back({network.S, -outputLiterals[0]});
				if (iff) {
					clauses.push_back({-network.S, outputLiterals[0]});
				}
			}
			else {
				assert(0);
			}
		}
	}
	else if (type == CONSTRAINT_ALL_LESS || type == CONSTRAINT_ALL_GREATER || type == CONSTRAINT_ALL_BOTH) {
		assert(0 < K && K == (int)outputLiterals.size());
		int k = 2;
		while (k < K) k *= 2;
		Formula F(freeVar);
		vector<int> input = inputLiterals;
		while ((int)input.size()%k != 0) {
			F.clauses.push_back({-F.freeVar});
			input.push_back({F.freeVar++});
		}
		bool d1 = (type == CONSTRAINT_ALL_LESS) || iff || (type == CONSTRAINT_ALL_BOTH);
		bool d2 = (type == CONSTRAINT_ALL_GREATER) || iff || (type == CONSTRAINT_ALL_BOTH);
		vector<int> output = Card(input, 0, (int)input.size()-1, k, F, d1, d2);
		clauses = F.clauses;
		for (int i = 0; i < K; i++) {
			if (d1) clauses.push_back({-output[i], outputLiterals[i]});
			if (d2) clauses.push_back({output[i], -outputLiterals[i]});
		}
	}
	else {
		assert(0);
	}
	for (auto& c : clauses) {
		for (int lit : addLiterals) {
			c.push_back(lit);
		}
	}
	return clauses;
}
}