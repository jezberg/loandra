#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>

#include "probleminstance.hpp"
#include "AMSLEX.hpp"
#include "timer.hpp"
#include "global.hpp"

#define F first
#define S second
#define getD(e,i) data[(e.B) + (i)]

using namespace std;
namespace maxPreprocessor {

AMSLEX::AMSLEX(const ProblemInstance& pi_) : pi(pi_) {
	data = 0;
	dataSize = 0;
}

AMSLEX::~AMSLEX() {
	free(data);
}

void AMSLEX::assumeSize(unsigned size) {
	if (dataSize < size) {
		unsigned newSize = max(2 * dataSize, size);
		data = (int*)realloc(data, newSize * sizeof(int));
		dataSize = newSize;
	}
}

bool AMSLEX::isPrefix(const vector<int>& a, const vector<int>& b) {
	for (unsigned i = 0; i < a.size(); i++) {
		if (a[i] != b[i]) return false;
	}
	return true;
}

bool AMSLEX::isPrefix(const vecP a, const vecP b) {
	for (unsigned i = 0; i < a.size(); i++) {
		if (getD(a, i) != getD(b, i)) return false;
	}
	return true;
}

bool AMSLEX::CSO1(const vector<vecP>& D, unsigned b, unsigned e, const vecP S, 
					   unsigned j, unsigned d, const vector<int>& d0Index) {
	while (j < S.size() && getD(S, j) < getD(D[b], d)) j++;
	if (j >= S.size()) return false;
	
	if (getD(S, j) == getD(D[b], d)) {
		unsigned ee = b;
		if (d == 0) {
			ee = d0Index[getD(S, j) + 1] - 1;
		}
		else if (e - b < BSconst) {
			while (ee + 1 <= e && getD(D[ee + 1], d) == getD(S, j)) ee++;
		}
		else {
			int mi = b;
			int ma = e - 1;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (getD(D[mid + 1], d) == getD(S, j)) {
					ee = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
		if (S.size() > d + 1 && D[b].size() == d + 1) {
			return true;
		}
		if (j + 1 <= S.size()) {
			if (CSO1(D, b, ee, S, j + 1, d + 1, d0Index)) {
				return true;
			}
		}
		b = ee + 1;
	}
	else {
		if (d == 0) {
			b = d0Index[getD(S, j)];
		}
		else if (e - b < BSconst) {
			while (b <= e && getD(D[b], d) < getD(S, j)) b++;
		}
		else {
			int mi = b;
			int ma = e;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (getD(D[mid], d) < getD(S, j)) {
					b = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
	}
	if (b <= e) return CSO1(D, b, e, S, j, d, d0Index);
	return false;
}

bool AMSLEX::CSO2(const vector<int>& D, unsigned b, unsigned e, const vector<int>& S, unsigned j, unsigned d) {
	while (j < S.size() && S[j] < pi.clauses[D[b]].lit[d]) j++;
	if (j >= S.size()) return false;
	
	if (S[j] == pi.clauses[D[b]].lit[d]) {
		unsigned ee = b;
		if (e - b < BSconst) {
			while (ee + 1 <= e && pi.clauses[D[ee + 1]].lit[d] == S[j]) ee++;
		}
		else {
			int mi = b;
			int ma = e - 1;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (pi.clauses[D[mid + 1]].lit[d] == S[j]) {
					ee = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
		if (S.size() > d + 1 && pi.clauses[D[b]].lit.size() == d + 1) {
			return true;
		}
		if (j + 1 <= S.size()) {
			if (CSO2(D, b, ee, S, j + 1, d + 1)) {
				return true;
			}
		}
		b = ee + 1;
	}
	else {
		if (e - b < BSconst) {
			while (b <= e && pi.clauses[D[b]].lit[d] < S[j]) b++;
		}
		else {
			int mi = b;
			int ma = e;
			while (mi <= ma) {
				int mid = (mi+ma)/2;
				if (pi.clauses[D[mid]].lit[d] < S[j]) {
					b = mid + 1;
					mi = mid + 1;
				}
				else {
					ma = mid - 1;
				}
			}
		}
	}
	if (b <= e) return CSO2(D, b, e, S, j, d);
	return false;
}

vector<int> AMSLEX::amsLexSENonPerm(const vector<int>& clauses) {
	vector<int> D(clauses.size());
	for (unsigned i = 0; i < clauses.size(); i++) {
		D[i] = clauses[i];
	}
	auto cmp = [&](int a, int b) {
		return pi.clauses[a].lit < pi.clauses[b].lit;
	};
	sort(D.begin(), D.end(), cmp);
	vector<int> subsumed(D.size());
	unsigned S = 0;
	for (unsigned i = 1; i < D.size(); i++) {
		if (pi.clauses[D[S]].lit.size() <= pi.clauses[D[i]].lit.size() && isPrefix(pi.clauses[D[S]].lit, pi.clauses[D[i]].lit)) {
			subsumed[i] = true;
		}
		else {
			S = i;
		}
	}
	for (unsigned i = 0; i + 1 < D.size(); i++) {
		if (!subsumed[i] && CSO2(D, i + 1, D.size() - 1, pi.clauses[D[i]].lit, 0, 0)) {
			subsumed[i] = true;
		}
	}
	vector<int> ret;
	for (unsigned i = 0; i < D.size(); i++) {
		if (subsumed[i]) ret.push_back(D[i]);
	}
	return ret;
}

vector<int> AMSLEX::amsLexSEPerm(const vector<int>& clauses) {
	vector<vecP> D(clauses.size());
	vector<int> d0Index;
	ALIt++;
	int pid = 0;
	int maxFreq = 0;
	unsigned dataP = 0;
	for (unsigned i = 0; i < clauses.size(); i++) {
		assumeSize(dataP + pi.clauses[clauses[i]].lit.size());
		memcpy(data + dataP, pi.clauses[clauses[i]].lit.data(), pi.clauses[clauses[i]].lit.size()*sizeof(int));
		D[i].B = dataP;
		D[i].E = dataP + pi.clauses[clauses[i]].lit.size();
		D[i].O = clauses[i];
		dataP += D[i].size();
		
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			if ((int)ALU.size() <= data[j]) {
				ALU.resize(data[j] + 1);
				ALI.resize(data[j] + 1);
			}
			if (ALU[data[j]] != ALIt) {
				ALU[data[j]] = ALIt;
				if (pid >= (int)ALF.size()) ALF.resize(pid + 1);
				ALF[pid] = 0;
				ALI[data[j]] = pid++;
			}
			data[j] = ALI[data[j]];
			ALF[data[j]]++;
			maxFreq = max(maxFreq, ALF[data[j]]);
		}
	}
	vector<vector<int> > cSort(maxFreq + 1);
	for (int i = 0; i < pid; i++) {
		cSort[ALF[i]].push_back(i);
	}
	vector<unsigned> perm(pid);
	unsigned i2 = 0;
	for (int i = 0; i <= maxFreq; i++) {
		for (int t : cSort[i]) {
			perm[t] = i2++;
		}
	}
	for (unsigned i = 0; i < D.size(); i++) {
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			data[j] = perm[data[j]];
		}
		sort(data + D[i].B, data + D[i].E);
	}
	auto cmp = [&](vecP a, vecP b) {
		if (getD(a, 0) < getD(b, 0)) return true;
		else if (getD(a, 0) > getD(b, 0)) return false;
		unsigned sz = min(a.size(), b.size());
		for (unsigned i = 1; i < sz; i++) {
			if (getD(a, i) < getD(b, i)) return true;
			else if (getD(a, i) > getD(b, i)) return false;
		}
		return a.size() < b.size();
	};
	vector<vector<vecP> > cSort2(pid);
	for (unsigned i = 0; i < D.size(); i++) {
		cSort2[getD(D[i], 0)].push_back(D[i]);
	}
	unsigned di = 0;
	for (int i = 0; i < pid; i++) {
		sort(cSort2[i].begin(), cSort2[i].end(), cmp);
		for (unsigned j = 0; j < cSort2[i].size(); j++) {
			D[di++] = cSort2[i][j];
		}
	}
	d0Index.resize(pid + 1);
	int idx = 0;
	for (unsigned i = 0; i < D.size(); i++) {
		if (i == 0 || getD(D[i], 0) > getD(D[i - 1], 0)) {
			while (idx <= getD(D[i], 0)) {
				d0Index[idx++] = i;
			}
		}
	}
	while (idx <= pid) {
		d0Index[idx++] = (int)D.size();
	}
	vector<int> subsumed(D.size());
	unsigned S = 0;
	for (unsigned i = 1; i < D.size(); i++) {
		if (D[S].size() <= D[i].size() && isPrefix(D[S], D[i])) {
			subsumed[i] = true;
		}
		else {
			S = i;
		}
	}
	for (unsigned i = 0; i + 1 < D.size(); i++) {
		if (!subsumed[i] && CSO1(D, i + 1, D.size() - 1, D[i], 0, 0, d0Index)) {
			subsumed[i] = true;
		}
	}
	vector<int> ret;
	for (unsigned i = 0; i < D.size(); i++) {
		if (subsumed[i]) ret.push_back(D[i].O);
	}
	return ret;
}

vector<int> AMSLEX::amsLexSE(const vector<int>& clauses) {
	if (clauses.size() <= 700) {
		return amsLexSENonPerm(clauses);
	}
	else {
		return amsLexSEPerm(clauses);
	}
}

pair<vector<int>, vector<int> > AMSLEX::amsLexSSRPerm(const vector<int>& c1, const vector<int>& c2, int var) {
	vector<vecP> D(c1.size() + c2.size());
	vector<int> d0Index;
	ALIt++;
	int pid = 0;
	int maxFreq = 0;
	unsigned dataP = 0;
	for (unsigned i = 0; i < c1.size(); i++) {
		if (pi.clauses[c1[i]].lit.size() == 1) {
			pair<vector<int>, vector<int> > ret;
			for (int a : c1) {
				ret.F.push_back(a);
			}
			for (int b : c2) {
				ret.S.push_back(b);
			}
			return ret;
		}
		assumeSize(dataP + pi.clauses[c1[i]].lit.size());
		memcpy(data + dataP, pi.clauses[c1[i]].lit.data(), pi.clauses[c1[i]].lit.size()*sizeof(int));
		D[i].B = dataP;
		D[i].E = dataP + pi.clauses[c1[i]].lit.size();
		D[i].O = 2*c1[i];
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			if (litVariable(data[j]) == var) {
				swap(data[j], data[D[i].E - 1]);
				break;
			}
		}
		D[i].E--;
		dataP += D[i].size();
		
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			if ((int)ALU.size() <= data[j]) {
				ALU.resize(data[j] + 1);
				ALI.resize(data[j] + 1);
			}
			if (ALU[data[j]] != ALIt) {
				ALU[data[j]] = ALIt;
				if (pid >= (int)ALF.size()) ALF.resize(pid + 1);
				ALF[pid] = 0;
				ALI[data[j]] = pid++;
			}
			data[j] = ALI[data[j]];
			ALF[data[j]]++;
			maxFreq = max(maxFreq, ALF[data[j]]);
		}
	}
	for (unsigned i = c1.size(); i < c1.size()+c2.size(); i++) {
		if (pi.clauses[c2[i - c1.size()]].lit.size() == 1) {
			pair<vector<int>, vector<int> > ret;
			for (int a : c1) {
				ret.F.push_back(a);
			}
			for (int b : c2) {
				ret.S.push_back(b);
			}
			return ret;
		}
		assumeSize(dataP + pi.clauses[c2[i - c1.size()]].lit.size());
		memcpy(data + dataP, pi.clauses[c2[i - c1.size()]].lit.data(), pi.clauses[c2[i - c1.size()]].lit.size()*sizeof(int));
		D[i].B = dataP;
		D[i].E = dataP + pi.clauses[c2[i - c1.size()]].lit.size();
		D[i].O = 2*c2[i - c1.size()] + 1;
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			if (litVariable(data[j]) == var) {
				swap(data[j], data[D[i].E - 1]);
				break;
			}
		}
		D[i].E--;
		dataP += D[i].size();
		
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			if ((int)ALU.size() <= data[j]) {
				ALU.resize(data[j] + 1);
				ALI.resize(data[j] + 1);
			}
			if (ALU[data[j]] != ALIt) {
				ALU[data[j]] = ALIt;
				if (pid >= (int)ALF.size()) ALF.resize(pid + 1);
				ALF[pid] = 0;
				ALI[data[j]] = pid++;
			}
			data[j] = ALI[data[j]];
			ALF[data[j]]++;
			maxFreq = max(maxFreq, ALF[data[j]]);
		}
	}
	vector<vector<int> > cSort(maxFreq + 1);
	for (int i = 0; i < pid; i++) {
		cSort[ALF[i]].push_back(i);
	}
	vector<unsigned> perm(pid);
	unsigned i2 = 0;
	for (int i = 0; i <= maxFreq; i++) {
		for (int t : cSort[i]) {
			perm[t] = i2++;
		}
	}
	for (unsigned i = 0; i < D.size(); i++) {
		for (unsigned j = D[i].B; j < D[i].E; j++) {
			data[j] = perm[data[j]];
		}
		sort(data + D[i].B, data + D[i].E);
	}
	auto cmp = [&](vecP a, vecP b) {
		if (getD(a, 0) < getD(b, 0)) return true;
		else if (getD(a, 0) > getD(b, 0)) return false;
		unsigned sz = min(a.size(), b.size());
		for (unsigned i = 1; i < sz; i++) {
			if (getD(a, i) < getD(b, i)) return true;
			else if (getD(a, i) > getD(b, i)) return false;
		}
		return a.size() < b.size();
	};
	vector<vector<vecP> > cSort2(pid);
	for (unsigned i = 0; i < D.size(); i++) {
		cSort2[getD(D[i], 0)].push_back(D[i]);
	}
	unsigned di = 0;
	for (int i = 0; i < pid; i++) {
		sort(cSort2[i].begin(), cSort2[i].end(), cmp);
		for (unsigned j = 0; j < cSort2[i].size(); j++) {
			D[di++] = cSort2[i][j];
		}
	}
	d0Index.resize(pid + 1);
	int idx = 0;
	for (unsigned i = 0; i < D.size(); i++) {
		if (i == 0 || getD(D[i], 0) > getD(D[i - 1], 0)) {
			while (idx <= getD(D[i], 0)) {
				d0Index[idx++] = i;
			}
		}
	}
	while (idx <= pid) {
		d0Index[idx++] = (int)D.size();
	}
	vector<int> subsumed(D.size());
	unsigned S = 0;
	for (unsigned i = 1; i < D.size(); i++) {
		if (D[S].size() <= D[i].size() && isPrefix(D[S], D[i])) {
			subsumed[i] = true;
		}
		else {
			S = i;
		}
	}
	for (unsigned i = 0; i + 1 < D.size(); i++) {
		if (!subsumed[i] && CSO1(D, i + 1, D.size() - 1, D[i], 0, 0, d0Index)) {
			subsumed[i] = true;
		}
	}
	pair<vector<int>, vector<int> > ret;
	for (unsigned i = 0; i < c1.size() + c2.size(); i++) {
		if (subsumed[i]) {
			if (D[i].O%2 == 0) {
				ret.F.push_back(D[i].O/2);
			}
			else {
				ret.S.push_back(D[i].O/2);
			}
		}
	}
	return ret;
}

pair<vector<int>, vector<int> > AMSLEX::amsLexSSR(const vector<int>& c1, const vector<int>& c2, int var) {
	return amsLexSSRPerm(c1, c2, var);
}

}