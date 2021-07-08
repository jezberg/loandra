vector<uint64_t> Preprocessor::getBVEHash(const vector<int>& cs, int var, int sw) const {
	vector<uint64_t> h(cs.size());
	for (unsigned i = 0; i < cs.size(); i++) {
		for (int l : pi.clauses[cs[i]].lit) {
			if (litVariable(l) != var) h[i] |= ((uint64_t)1<<(uint64_t)((l^sw)&63));
		}
	}
	return h;
}

pair<vector<int>, int> Preprocessor::searchAndOr(int dLit, int defVars, const set<int>& binaryClauseLits) const {
	vector<int> defClauses;
	for (int c : pi.litClauses[dLit]) {
		if (pi.clauses[c].lit.size() > 2 && ((int)pi.clauses[c].lit.size() - 1 < defVars || defVars == 0)) {
			bool f = false;
			for (int lit : pi.clauses[c].lit) {
				if (lit != dLit) {
					if (binaryClauseLits.count(litNegation(lit)) == 0) {
						f = true;
						break;
					}
				}
			}
			if (!f) {
				defVars = pi.clauses[c].lit.size() - 1;
				defClauses.clear();
				defClauses.push_back(c);
				for (int cc : pi.litClauses[litNegation(dLit)]) {
					if (pi.clauses[cc].lit.size() == 2) {
						int tl;
						if (litNegation(pi.clauses[cc].lit[0]) == dLit) {
							tl = pi.clauses[cc].lit[1];
						}
						else {
							tl = pi.clauses[cc].lit[0];
						}
						for (int lit : pi.clauses[c].lit) {
							if (litNegation(lit) == tl) {
								defClauses.push_back(cc);
								break;
							}
						}
					}
				}
				assert(defClauses.size() >= pi.clauses[c].lit.size());
			}
		}
	}
	return {defClauses, defVars};
}

pair<vector<int>, int> Preprocessor::searchITE(int var) const {
	vector<int> defClauses;
	int p3 = 0;
	int n3 = 0;
	for (int c : pi.litClauses[posLit(var)]) {
		if (pi.clauses[c].lit.size() == 3) {
			p3++;
			if (p3 >= 2) break;
		}
	}
	for (int c : pi.litClauses[negLit(var)]) {
		if (pi.clauses[c].lit.size() == 3) {
			n3++;
			if (n3 >= 2) break;
		}
	}
	if (p3 < 2 || n3 < 2) return {vector<int>(), 0};
	map<int, map<int, int> > pMap;
	map<int, map<int, int> > nMap;
	for (int c1 : pi.litClauses[posLit(var)]) {
		if (pi.clauses[c1].lit.size() == 3) {
			int l1 = pi.clauses[c1].lit[0];
			int l2 = pi.clauses[c1].lit[1];
			int l3 = pi.clauses[c1].lit[2];
			if (litVariable(l1) == var) swap(l1, l3);
			if (litVariable(l2) == var) swap(l2, l3);
			assert(litVariable(l1) != var && litVariable(l2) != var);
			pMap[l1][l2] = c1;
			pMap[l2][l1] = c1;
		}
	}
	for (int c1 : pi.litClauses[negLit(var)]) {
		if (pi.clauses[c1].lit.size() == 3) {
			int l1 = pi.clauses[c1].lit[0];
			int l2 = pi.clauses[c1].lit[1];
			int l3 = pi.clauses[c1].lit[2];
			if (litVariable(l1) == var) swap(l1, l3);
			if (litVariable(l2) == var) swap(l2, l3);
			assert(litVariable(l1) != var && litVariable(l2) != var);
			nMap[l1][l2] = c1;
			nMap[l2][l1] = c1;
		}
	}
	for (auto k : pMap) {
		int c = k.F;
		pair<int, int> tf = {-1, -1};
		pair<int, int> ff = {-1, -1};
		for (auto h : pMap[c]) {
			if (pMap[c].count(h.F) > 0 && nMap[c].count(litNegation(h.F)) > 0) {
				tf = {pMap[c][h.F], nMap[c][litNegation(h.F)]};
			}
		}
		if (tf.F == -1) continue;
		for (auto h : pMap[litNegation(c)]) {
			if (pMap[litNegation(c)].count(h.F) > 0 && nMap[litNegation(c)].count(litNegation(h.F)) > 0) {
				ff = {pMap[litNegation(c)][h.F], nMap[litNegation(c)][litNegation(h.F)]};
			}
		}
		if (tf.F >= 0 && ff.F >= 0) {
			return {{tf.F, tf.S, ff.F, ff.S}, 4};
		}
	}
	return {vector<int>(), 0};
}

pair<vector<int>, int> Preprocessor::searchXor(int var) const {
	if (pi.litClauses[posLit(var)].size() < 2 || pi.litClauses[negLit(var)].size() < 2) return {vector<int>(), 0};
	int maxXorVar = 2;
	while ((1 << maxXorVar) <= (int)min(pi.litClauses[posLit(var)].size(), pi.litClauses[negLit(var)].size())) maxXorVar++;
	vector<unordered_map<int, int> > pCnt(maxXorVar + 2);
	vector<unordered_map<int, int> > nCnt(maxXorVar + 2);
	
	for (int c : pi.litClauses[posLit(var)]) {
		if ((int)pi.clauses[c].lit.size() > maxXorVar + 1) continue;
		if (pi.clauses[c].lit.size() < 3) continue;
		for (int l : pi.clauses[c].lit) {
			pCnt[(int)pi.clauses[c].lit.size() - 1][l]++;
		}
	}
	for (int c : pi.litClauses[negLit(var)]) {
		if ((int)pi.clauses[c].lit.size() > maxXorVar + 1) continue;
		if (pi.clauses[c].lit.size() < 3) continue;
		for (int l : pi.clauses[c].lit) {
			nCnt[(int)pi.clauses[c].lit.size() - 1][l]++;
		}
	}
	vector<pair<vector<int>, int> > cand;
	for (int c : pi.litClauses[posLit(var)]) {
		if ((int)pi.clauses[c].lit.size() > maxXorVar + 1) continue;
		if (pi.clauses[c].lit.size() < 3) continue;
		int sz = pi.clauses[c].lit.size() - 1;
		bool f = false;
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) == var) continue;
			if (pCnt[sz][l] < (1 << (sz-2))) {
				f = true;
				break;
			}
			if (nCnt[sz][l] < (1 << (sz-2))) {
				f = true;
				break;
			}
		}
		if (f) continue;
		cand.push_back({vector<int>(), c});
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) == var) continue;
			cand.back().F.push_back(litVariable(l));
		}
	}
	for (int c : pi.litClauses[negLit(var)]) {
		if ((int)pi.clauses[c].lit.size() > maxXorVar + 1) continue;
		if (pi.clauses[c].lit.size() < 3) continue;
		int sz = pi.clauses[c].lit.size() - 1;
		bool f = false;
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) == var) continue;
			if (pCnt[sz][l] < (1 << (sz-2))) {
				f = true;
				break;
			}
			if (nCnt[sz][l] < (1 << (sz-2))) {
				f = true;
				break;
			}
		}
		if (f) continue;
		cand.push_back({vector<int>(), c});
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) == var) continue;
			cand.back().F.push_back(litVariable(l));
		}
	}
	sort(cand.begin(), cand.end());
	vector<int> tCls;
	for (unsigned i = 0; i < cand.size(); i++) {
		if (i == 0 || cand[i].F != cand[i - 1].F) {
			tCls.clear();
		}
		tCls.push_back(cand[i].S);
		if (i + 1 == cand.size() || cand[i + 1].F != cand[i].F) {
			int sz = cand[i].F.size();
			if ((1 << sz) > (int)tCls.size()) continue;
			vector<int> enC;
			vector<int> onC;
			vector<int> epC;
			vector<int> opC;
			for (int b = 0; b < (1<<(int)cand[i].F.size()); b++) {
				bool ff = false;
				int enf = -1;
				int onf = -1;
				int epf = -1;
				int opf = -1;
				for (int c : tCls) {
					bool f = false;
					bool pol = false;
					int pt = 0;
					for (int l : pi.clauses[c].lit) {
						if (litValue(l) == false) pt ^= 1;
						if (l == posLit(var)) pol = true;
						else if (l == negLit(var)) pol = false;
						else {
							bool lf = false;
							int lv = litVariable(l);
							for (int j = 0; j < (int)cand[i].F.size(); j++) {
								if (lv == cand[i].F[j]) {
									if (litValue(l) == false && (b&(1<<j)) == 0) {
										lf = true;
										break;
									}
									if (litValue(l) == true && (b&(1<<j)) > 0) {
										lf = true;
										break;
									}
								}
							}
							if (!lf) {
								f = true;
								break;
							}
						}
					}
					if (!f) {
						ff = true;
						if (pol == false && pt == 0) enf = c;
						else if (pol == false && pt == 1) onf = c;
						else if (pol == true && pt == 0) epf = c;
						else if (pol == true && pt == 1) opf = c;
					}
				}
				if (enf >= 0) {
					enC.push_back(enf);
					assert(epf == -1);
				}
				if (onf >= 0) {
					onC.push_back(onf);
					assert(opf == -1);
				}
				if (epf >= 0) {
					epC.push_back(epf);
					assert(enf == -1);
				}
				if (opf >= 0) {
					opC.push_back(opf);
					assert(onf == -1);
				}
				if (!ff) break;
			}
			if ((int)enC.size() == (1 << sz)/2 && (int)epC.size() == (1 << sz)/2) {
				vector<int> ret;
				for (int c : enC) {
					ret.push_back(c);
				}
				for (int c : epC) {
					ret.push_back(c);
				}
				assert((int)ret.size() == (1 << sz));
				return {ret, (int)ret.size()};
			}
			if ((int)onC.size() == (1 << sz)/2 && (int)opC.size() == (1 << sz)/2) {
				vector<int> ret;
				for (int c : onC) {
					ret.push_back(c);
				}
				for (int c : opC) {
					ret.push_back(c);
				}
				assert((int)ret.size() == (1 << sz));
				return {ret, (int)ret.size()};
			}
		}
	}
	return {vector<int>(), 0};
}

vector<int> Preprocessor::tryBVEGE(int var) {
	set<int> negBinaryClauseLits, posBinaryClauseLits;
	
	for (int c : pi.litClauses[negLit(var)]) {
		if (pi.clauses[c].lit.size() == 2) {
			if (litVariable(pi.clauses[c].lit[0]) == var) {
				negBinaryClauseLits.insert(pi.clauses[c].lit[1]);
			}
			else {
				negBinaryClauseLits.insert(pi.clauses[c].lit[0]);
			}
		}
	}
	for (int c : pi.litClauses[posLit(var)]) {
		if (pi.clauses[c].lit.size() == 2) {
			if (litVariable(pi.clauses[c].lit[0]) == var) {
				posBinaryClauseLits.insert(pi.clauses[c].lit[1]);
			}
			else {
				posBinaryClauseLits.insert(pi.clauses[c].lit[0]);
			}
		}
	}
	
	//Find the shortest definition
	vector<int> defClauses;
	int defVars = 0;
	if (defVars == 0) {
		auto tr = searchXor(var);
		if (tr.S > 0) {
			defVars = tr.S;
			defClauses = tr.F;
		}
	}
	if (defVars == 0) {
		auto tr = searchAndOr(posLit(var), defVars, negBinaryClauseLits);
		if (tr.S > 0) {
			defVars = tr.S;
			defClauses = tr.F;
		}
	}
	if (defVars == 0) {
		auto tr = searchAndOr(negLit(var), defVars, posBinaryClauseLits);
		if (tr.S > 0) {
			defVars = tr.S;
			defClauses = tr.F;
		}
	}
	if (defVars == 0) {
		auto tr = searchITE(var);
		if (tr.S > 0) {
			defVars = tr.S;
			defClauses = tr.F;
		}
	}
	if (defVars > 0) {
		return defClauses;
	}
	return vector<int>();
}

// Dont give labels to this
int Preprocessor::tryBVE2(int var) {
	int sizeLimit = pi.litClauses[posLit(var)].size() + pi.litClauses[negLit(var)].size();
	
	vector<int> isPosDef, isNegDef;
	vector<int> defClauses;
	if (BVEgate) {
		defClauses = tryBVEGE(var);
	}
	bool defFound = false;
	if (defClauses.size() > 0) {
		defFound = true;
		rLog.gatesExtracted++;
		isPosDef.resize(pi.litClauses[posLit(var)].size());
		isNegDef.resize(pi.litClauses[negLit(var)].size());
		sort(defClauses.begin(), defClauses.end());
		for (int i = 0; i < (int)pi.litClauses[posLit(var)].size(); i++) {
			if (binary_search(defClauses.begin(), defClauses.end(), pi.litClauses[posLit(var)][i])) {
				isPosDef[i] = 1;
			}
		}
		for (int i = 0; i < (int)pi.litClauses[negLit(var)].size(); i++) {
			if (binary_search(defClauses.begin(), defClauses.end(), pi.litClauses[negLit(var)][i])) {
				isNegDef[i] = 1;
			}
		}
	}
	
	//the condition sizelimit >= 6 doesnt really help
	if (sizeLimit >= 6) { // magic constant
		vector<uint64_t> h1=getBVEHash(pi.litClauses[posLit(var)], var, 0);
		vector<uint64_t> h2=getBVEHash(pi.litClauses[negLit(var)], var, 1);
		int hcnt = 0;
		for (int i = 0; i < (int)pi.litClauses[posLit(var)].size(); i++) {
			for (int ii = 0; ii < (int)pi.litClauses[negLit(var)].size(); ii++) {
				if (defFound && (isPosDef[i] == isNegDef[ii])) continue;
				if ((h1[i] & h2[ii]) == 0) {
					hcnt++;
					if (hcnt > sizeLimit) return 0;
				}
			}
		}
	}
	
	vector<vector<int> > newClauses;
	for (int i = 0; i < (int)pi.litClauses[posLit(var)].size(); i++) {
		for (int ii = 0; ii < (int)pi.litClauses[negLit(var)].size(); ii++) {
			if (defFound && (isPosDef[i] == isNegDef[ii])) continue;
			const vector<int>& c1 = pi.clauses[pi.litClauses[posLit(var)][i]].lit;
			const vector<int>& c2 = pi.clauses[pi.litClauses[negLit(var)][ii]].lit;
			unsigned j2 = 0;
			bool f = false;
			for (unsigned j = 0; j < c1.size(); j++) {
				while (j2 < c2.size() && c2[j2] <= c1[j]) {
					if (litNegation(c2[j2]) == c1[j] && litVariable(c1[j]) != var) {
						f = true;
						break;
					}
					j2++;
				}
				if (f) break;
				if (j2 < c2.size() && litNegation(c2[j2]) == c1[j] && litVariable(c1[j]) != var) {
					f = true;
					break;
				}
			}
			
			if (f) continue;
			if ((int)newClauses.size() >= sizeLimit) return 0;
			newClauses.push_back(vector<int>());
			j2 = 0;
			for (unsigned j = 0; j < c1.size(); j++) {
				while (j2 < c2.size() && c2[j2] <= c1[j]) {
					if (litVariable(c2[j2]) != var) {
						newClauses.back().push_back(c2[j2]);
					}
					j2++;
				}
				if (!(newClauses.back().size() > 0 && newClauses.back().back() == c1[j]) && litVariable(c1[j]) != var) {
					newClauses.back().push_back(c1[j]);
				}
			}
			while (j2 < c2.size()) {
				if (litVariable(c2[j2]) != var) {
					newClauses.back().push_back(c2[j2]);
				}
				j2++;
			}
		}
	}
	
	vector<vector<int> > nClauses;
	vector<int> toRemove;
	for (int c : pi.litClauses[posLit(var)]) {
		toRemove.push_back(c);
	}
	for (int c : pi.litClauses[negLit(var)]) {
		toRemove.push_back(c);
		nClauses.push_back(pi.clauses[c].lit);
	}
	for (int c : toRemove) {
		pi.removeClause(c);
	}
	for (auto& nc : newClauses) {
		pi.addClause(nc);
	}
	trace.BVE(var, nClauses);
	
	rLog.removeVariable(1);
	rLog.removeClause((int)toRemove.size() - (int)newClauses.size());
	
	assert(pi.litClauses[posLit(var)].size() == 0 && pi.litClauses[negLit(var)].size() == 0);
	return 1;
}

// Dont give labels to this
int Preprocessor::tryBVE(int var) {
	vector<vector<int> > newClauses;
	int sizeLimit = pi.litClauses[posLit(var)].size() + pi.litClauses[negLit(var)].size();
	for (int i = 0; i < (int)pi.litClauses[posLit(var)].size(); i++) {
		for (int ii = 0; ii < (int)pi.litClauses[negLit(var)].size(); ii++) {
			// Check if resulting clause will be a tautology
			Clause& c1 = pi.clauses[pi.litClauses[posLit(var)][i]];
			Clause& c2 = pi.clauses[pi.litClauses[negLit(var)][ii]];
			int j2 = 0;
			bool f = false;
			for (int j = 0; j < (int)c1.lit.size(); j++) {
				while (j2 < (int)c2.lit.size() && c2.lit[j2] <= c1.lit[j]) {
					if (litNegation(c2.lit[j2]) == c1.lit[j] && litVariable(c1.lit[j]) != var) {
						f = true;
						break;
					}
					j2++;
				}
				if (f) break;
				if (j2 < (int)c2.lit.size() && litNegation(c2.lit[j2]) == c1.lit[j] && litVariable(c1.lit[j]) != var) {
					f = true;
					break;
				}
			}
			if (f) continue;
			if ((int)newClauses.size() >= sizeLimit) return 0;
			newClauses.push_back(vector<int>());
			j2 = 0;
			for (int j = 0; j < (int)c1.lit.size(); j++) {
				while ((int)j2 < (int)c2.lit.size() && c2.lit[j2] <= c1.lit[j]) {
					if (litVariable(c2.lit[j2]) != var) {
						newClauses.back().push_back(c2.lit[j2]);
					}
					j2++;
				}
				if (!((int)newClauses.back().size() > 0 && newClauses.back().back() == c1.lit[j]) && litVariable(c1.lit[j]) != var) {
					newClauses.back().push_back(c1.lit[j]);
				}
			}
			while ((int)j2 < (int)c2.lit.size()) {
				if (litVariable(c2.lit[j2]) != var) {
					newClauses.back().push_back(c2.lit[j2]);
				}
				j2++;
			}
			
			// Check if valid, this can be deleted
			for (int j = 1; j < (int)newClauses.back().size(); j++) {
				assert(litVariable(newClauses.back()[j]) != var);
				assert(litVariable(newClauses.back()[j - 1]) != var);
				assert(newClauses.back()[j] > newClauses.back()[j - 1]);
				assert(litNegation(newClauses.back()[j]) != newClauses.back()[j - 1]);
			}
		}
	}
	
	vector<vector<int> > nClauses;
	
	vector<int> toRemove;
	for (int c : pi.litClauses[posLit(var)]) {
		toRemove.push_back(c);
	}
	for (int c : pi.litClauses[negLit(var)]) {
		toRemove.push_back(c);
		nClauses.push_back(pi.clauses[c].lit);
	}
	for (int c : toRemove) {
		pi.removeClause(c);
	}
	for (int i = 0; i < (int)newClauses.size(); i++) {
		pi.addClause(newClauses[i]);
	}
	trace.BVE(var, nClauses);
	
	rLog.removeVariable(1);
	rLog.removeClause((int)toRemove.size() - (int)newClauses.size());
	
	//Check if valid, this can be deleted
	assert(pi.litClauses[posLit(var)].size() == 0 && pi.litClauses[negLit(var)].size() == 0);
	return 1;
}

int Preprocessor::doBVE() {
	rLog.startTechnique(Log::Technique::BVE);
	if (!rLog.requestTime(Log::Technique::BVE)) {
		rLog.stopTechnique(Log::Technique::BVE);
		return 0;
	}
	int eliminated = 0;
	vector<int> checkVar = pi.tl.getTouchedVariables("BVE");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[negLit(var1)].size() + pi.litClauses[posLit(var1)].size() < pi.litClauses[negLit(var2)].size() + pi.litClauses[posLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	bool skip = false;
	if (skipTechnique > 0 && (int)checkVar.size() >= skipTechnique * 4) {
		for (int tc = 0; tc < skipTechnique; tc++) {
			if (!rLog.requestTime(Log::Technique::BVE)) break;
			int var = checkVar[getRand(0, (int)checkVar.size() - 1)];
			if (pi.isLabel[var] == 0 && (pi.litClauses[posLit(var)].size() > 0 || pi.litClauses[negLit(var)].size() > 0)) {
				eliminated += tryBVE2(var);
			}
		}
		if (eliminated == 0) {
			skip = true;
			log("BVE skipped");
		}
	}
	if (!skip) {
		for (int var : checkVar) {
			if (!rLog.requestTime(Log::Technique::BVE)) break;
			if (pi.isLabel[var] == 0 && (pi.litClauses[posLit(var)].size() > 0 || pi.litClauses[negLit(var)].size() > 0)) {
				eliminated += tryBVE2(var);
			}
		}
	}
	log(eliminated, " eliminated by BVE");
	rLog.stopTechnique(Log::Technique::BVE);
	return eliminated;
}

void Preprocessor::doBVE2() {
	for (int var = 0; var < pi.vars; var++) {
		if (pi.isLabel[var] == 0 && (pi.litClauses[posLit(var)].size() > 0 || pi.litClauses[negLit(var)].size() > 0)) {
			if (tryBVE(var) != 0) {
				print("fail BVE ", var + 1);
				abort();
			}
		}
	}
}