void Preprocessor::addBVAHash(vector<int>& lits, unordered_map<uint64_t, int>& hashes) {
	if (lits.size() <= 1) return;
	if (sfH.size() < lits.size() + 1) sfH.resize(lits.size() + 1);
	if (tMul.size() < lits.size() + 1) tMul.resize(lits.size() + 1);
	tMul[lits.size()] = 1;
	sfH[lits.size()] = 0;
	uint64_t mul = 1;
	for (int i = (int)lits.size() - 1; i >= 0; i--) {
		sfH[i] = sfH[i + 1];
		sfH[i] += mul * (uint64_t)lits[i];
		mul *= polyHashMul;
		tMul[i] = mul;
	}
	hashes[sfH[1]]++;
	uint64_t preH = 0;
	for (int i = 0; i < (int)lits.size() - 1; i++) {
		preH *= polyHashMul;
		preH += (uint64_t)lits[i];
		uint64_t mHash = preH * tMul[i + 2] + sfH[i + 2];
		hashes[mHash]++;
	}
}

int Preprocessor::canBVA(int c, int d, int lit) {
	if (pi.clauses[d].lit.size() != pi.clauses[c].lit.size()) return -1;
	int hc = -1;
	int hd = -1;
	unsigned i1 = 0;
	unsigned i2 = 0;
	while (i1 < pi.clauses[c].lit.size() || i2 < pi.clauses[d].lit.size()) {
		if (i1 < pi.clauses[c].lit.size() && i2 < pi.clauses[d].lit.size() && pi.clauses[c].lit[i1] == pi.clauses[d].lit[i2]) {
			i1++;
			i2++;
		}
		else if (i1 < pi.clauses[c].lit.size() && (i2 == pi.clauses[d].lit.size() || pi.clauses[d].lit[i2] > pi.clauses[c].lit[i1])) {
			if (hc != -1) {
				return -1;
			}
			hc = pi.clauses[c].lit[i1];
			if (hc != lit) {
				return -1;
			}
			i1++;
		}
		else if (i2 < pi.clauses[d].lit.size() && (i1 == pi.clauses[c].lit.size() || pi.clauses[c].lit[i1] > pi.clauses[d].lit[i2])) {
			if (hd != -1) {
				return -1;
			}
			hd = pi.clauses[d].lit[i2];
			i2++;
		}
	}
	return hd;
}

int Preprocessor::tryBVA(int lit, unordered_map<uint64_t, int>& hashes) {
	const bool noLabels = false;
	if (noLabels && pi.isLabel[litVariable(lit)]) return 0;
	removeDuplicateClauses();
	vector<int> mLit = {lit};
	vector<int> mCls = pi.litClauses[lit];
	int redu = 0;
	int taut = 0;
	while (1) {
		vector<pair<int, int> > P;
		for (int c : mCls) {
			if (pi.clauses[c].lit.size() == 1) continue;
			int lMin = pi.clauses[c].lit[0];
			if (lMin == lit) lMin = pi.clauses[c].lit[1];
			if (hashes.size() > 0) {
				uint64_t hash = 0;
				for (int l : pi.clauses[c].lit) {
					if (l != lit) {
						hash *= polyHashMul;
						hash += (uint64_t)l;
					}
				}
				int asd = hashes[hash];
				assert(asd > 0);
				if (asd == 1) continue;
			}
			for (unsigned i = 1; i < pi.clauses[c].lit.size(); i++) {
				if (pi.litClauses[pi.clauses[c].lit[i]].size() < pi.litClauses[lMin].size() && pi.clauses[c].lit[i] != lit) {
					lMin = pi.clauses[c].lit[i];
				}
			}
			for (int d : pi.litClauses[lMin]) {
				int hd = canBVA(c, d, lit);
				if (hd != -1) {
					assert(hd != lit);
					if (!(pi.isLabel[litVariable(hd)] && noLabels)) P.push_back({hd, c});
				}
			}
		}
		if (P.size() == 0) break;
		sort(P.begin(), P.end());
		P.erase(unique(P.begin(), P.end()), P.end());
		for (auto p : P) {
			if (p.F == litNegation(lit)) {
				int rClause = -1;
				for (int d : pi.litClauses[p.F]) {
					int hd = canBVA(p.S, d, lit);
					if (hd == p.F) {
						rClause = d;
						break;
					}
				}
				assert(rClause != -1);
				pi.removeClause(rClause);
				pi.removeLiteralFromClause(lit, p.S);
				if (hashes.size() > 0) addBVAHash(pi.clauses[p.S].lit, hashes);
				return 1;
			}
		}
		sort(mLit.begin(), mLit.end());
		int fq = 0;
		int lMax = -1;
		int fMax = 0;
		unsigned i2 = 0;
		for (unsigned i = 0; i < P.size(); i++) {
			if (i > 0 && P[i].F != P[i - 1].F) {
				fq = 0;
			}
			fq++;
			while (i2 < mLit.size() && mLit[i2] < P[i].F) {
				i2++;
			}
			if (i2 == mLit.size() || mLit[i2] != P[i].F) {
				if (fq > fMax) {
					fMax = fq;
					lMax = P[i].F;
				}
			}
		}
		if (lMax == -1) break;
		int nRedu = ((int)mLit.size() + 1) * fMax - ((int)mLit.size() + 1) - fMax - taut;
		if (nRedu <= redu) {
			// Try allowing tautologies in resolvents to find more
			fq = 0;
			int tMax = 0;
			fMax = -1;
			lMax = -1;
			int ft = 0;
			i2 = 0;
			for (unsigned i = 0; i < P.size(); i++) {
				if (i == 0 || P[i].F != P[i - 1].F) {
					fq = 0;
					ft = 0;
					for (int c : mCls) {
						if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litNegation(P[i].F))) {
							fq++;
							ft++;
						}
					}
				}
				fq++;
				while (i2 < mLit.size() && mLit[i2] < P[i].F) {
					i2++;
				}
				if (i2 == mLit.size() || mLit[i2] != P[i].F) {
					int nR = ((int)mLit.size() + 1) * fq - ((int)mLit.size() + 1) - fq - taut - ft;
					if (nR > nRedu) {
						nRedu = nR;
						lMax = P[i].F;
						fMax = fq;
						tMax = ft;
					}
				}
			}
			if (nRedu > redu) {
				redu = nRedu;
				mLit.push_back(lMax);
				vector<int> tmCls;
				for (int c : mCls) {
					if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litNegation(lMax))) {
						tmCls.push_back(c);
					}
				}
				mCls = tmCls;
				taut += tMax;
				for (unsigned i = 0; i < P.size(); i++) {
					if (P[i].F == lMax) {
						mCls.push_back(P[i].S);
					}
				}
				assert((int)mCls.size() == fMax);
			}
			else {
				break;
			}
		}
		else {
			redu = nRedu;
			mLit.push_back(lMax);
			mCls.clear();
			for (unsigned i = 0; i < P.size(); i++) {
				if (P[i].F == lMax) {
					mCls.push_back(P[i].S);
				}
			}
			assert((int)mCls.size() == fMax);
		}
	}
	if (redu == 0) return 0;
	assert(redu == (int)mLit.size() * (int)mCls.size() - (int)mLit.size() - (int)mCls.size() - taut);
	sort(mLit.begin(), mLit.end());
	vector<int> rmClauses;
	for (int c : mCls) {
		rmClauses.push_back(c);
		assert(pi.clauses[c].lit.size() > 1);
		int lMin = pi.clauses[c].lit[0];
		if (lMin == lit) lMin = pi.clauses[c].lit[1];
		for (unsigned i = 1; i < pi.clauses[c].lit.size(); i++) {
			if (pi.litClauses[pi.clauses[c].lit[i]].size() < pi.litClauses[lMin].size() && pi.clauses[c].lit[i] != lit) {
				lMin = pi.clauses[c].lit[i];
			}
		}
		for (int d : pi.litClauses[lMin]) {
			int hd = canBVA(c, d, lit);
			if (hd != -1 && binary_search(mLit.begin(), mLit.end(), hd)) {
				rmClauses.push_back(d);
			}
		}
	}
	sort(rmClauses.begin(), rmClauses.end());
	rmClauses.erase(unique(rmClauses.begin(), rmClauses.end()), rmClauses.end());
	assert((int)rmClauses.size() >= (int)mLit.size() * (int)mCls.size() - taut);
	assert((int)rmClauses.size() <= (int)mLit.size() * (int)mCls.size());
	int nVar = pi.addVar();
	int realRedu = 0;
	for (int l : mLit) {
		vector<int> newClause = {l, posLit(nVar)};
		pi.addClause(newClause);
		realRedu--;
		if (hashes.size() > 0) addBVAHash(newClause, hashes);
	}
	for (int c : mCls) {
		vector<int> newClause;
		for (int l : pi.clauses[c].lit) {
			if (!binary_search(mLit.begin(), mLit.end(), l)) {
				newClause.push_back(l);
			}
		}
		newClause.push_back(negLit(nVar));
		pi.addClause(newClause);
		realRedu--;
		if (hashes.size() > 0) addBVAHash(newClause, hashes);
	}
	for (int c : rmClauses) {
		pi.removeClause(c);
		realRedu++;
	}
	assert(realRedu >= redu);
	assert(realRedu <= redu + taut);
	rLog.removeClause(realRedu);
	rLog.removeVariable(-1);
	return realRedu;
}

int Preprocessor::doBVA() {
	rLog.startTechnique(Log::Technique::BVA);
	if (!rLog.requestTime(Log::Technique::BVA)) {
		rLog.stopTechnique(Log::Technique::BVA);
		return 0;
	}
	int removed = 0;
	vector<int> checkLit = pi.tl.getTouchedLiterals("BVA");
	vector<int> hClauses = pi.tl.getModClauses("BVAhash");
	for (int c : hClauses) {
		if (!rLog.requestTime(Log::Technique::BVA)) break;
		if (!pi.isClauseRemoved(c)) {
			addBVAHash(pi.clauses[c].lit, BVAHashTable);
		}
	}
	if (rLog.isTimeLimit()) {
		// The best literals to consider should be the ones that occur most
		// (at least the paper had this heuristic)
		auto cmp = [&](int lit1, int lit2) {
			return pi.litClauses[lit1].size() > pi.litClauses[lit2].size();
		};
		sort(checkLit.begin(), checkLit.end(), cmp);
	}
	for (int lit : checkLit) {
		if (!rLog.requestTime(Log::Technique::BVA)) break;
		if (pi.litClauses[lit].size() < 2) continue;
		removed += tryBVA(lit, BVAHashTable);
	}
	pi.tl.setItr("BVAhash");
	log(removed, " clauses removed by BVA");
	rLog.stopTechnique(Log::Technique::BVA);
	return removed;
}

void Preprocessor::doBVA2() {
	rLog.startTechnique(Log::Technique::BVA);
	unordered_map<uint64_t, int> hashes;
	for (int lit = 0; lit < 2*pi.vars; lit++) {
		if (tryBVA(lit, hashes)) {
			if (pi.litClauses[lit].size() < 2) continue;
			print("fail BVA ", litToDimacs(lit));
			abort();
		}
	}
	rLog.stopTechnique(Log::Technique::BVA);
}