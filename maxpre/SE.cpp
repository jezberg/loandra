// Is clause a subsumed by clause b?
// Supposes that a and b are sorted
bool Preprocessor::isSubsumed(const vector<int>& a, const vector<int>& b) const {
	unsigned i2 = 0;
	for (unsigned i = 0; i < b.size(); i++) {
		while (i2 < a.size() && a[i2] < b[i]) {
			i2++;
		}
		if (i2 >= a.size() || a[i2] != b[i]) return false;
	}
	return true;
}

// Slow implementation
int Preprocessor::trySESlow(int lit) {
	vector<int> toRemove;
	
	for (int c1 : pi.litClauses[lit]) {
		for (int c2 : pi.litClauses[lit]) {
			if (c1 <= c2) continue;
			if (pi.clauses[c1].lit.size() > pi.clauses[c2].lit.size() && pi.clauses[c2].isHard()) {
				if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2].lit)) {
					toRemove.push_back(c1);
				}
			}
			if (pi.clauses[c2].lit.size() > pi.clauses[c1].lit.size() && pi.clauses[c1].isHard()) {
				if (isSubsumed(pi.clauses[c2].lit, pi.clauses[c1].lit)) {
					toRemove.push_back(c2);
				}
			}
		}
	}
	int removed = 0;
	for (int c : toRemove) {
		if (!pi.isClauseRemoved(c)) {
			pi.removeClause(c);
			removed++;
		}
	}
	return removed;
}

void Preprocessor::trySEHash(vector<int>& clauses, int tLit, vector<int>& toRemove) {
	int k = 1;
	int n = clauses.size();
	while ((1<<k) < n) k++;
	vector<vector<pair<int, uint64_t> > > has(1 << k);
	for (int c : clauses) {
		uint64_t h = 0;
		int ml = pi.clauses[c].lit[0];
		for (int l : pi.clauses[c].lit) {
			h |= ((uint64_t)1 << (uint64_t)(l%k));
			if (pi.litClauses[l].size() < pi.litClauses[ml].size()) ml = l;
		}
		if (tLit == -1 || ml == tLit) {
			has[h].push_back({c, pi.clauses[c].hash});
		}
	}
	
	for (int c1 : clauses) {
		int64_t h = 0;
		uint64_t h64 = pi.clauses[c1].hash;
		for (int l : pi.clauses[c1].lit) {
			h |= ((int64_t)1 << (int64_t)(l%k));
		}
		for (int64_t sub = 0; (sub = (sub - h) & h);) {
			for (auto c2 : has[sub]) {
				if (!pi.canSubsume2(h64, c2.S)) continue;
				if (c1 == c2.F) continue;
				if (pi.clauses[c1].lit.size() == pi.clauses[c2.F].lit.size() && pi.canSubsume(c1, c2.F) && pi.canSubsume(c2.F, c1) && c1 > c2.F) {
					if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2.F].lit)) {
						if (pi.clauses[c1].isHard() && pi.clauses[c2.F].isHard()) {
							toRemove.push_back(c2.F);
						}
						else if(pi.clauses[c1].isHard()) {
							toRemove.push_back(c2.F);
						}
						else if(pi.clauses[c2.F].isHard()) {
							toRemove.push_back(c1);
						}
						else {
							pi.clauses[c1].weight += pi.clauses[c2.F].weight;
							pi.clauses[c2.F].weight = 0;
							toRemove.push_back(c2.F);
						}
					}
				}
				else if (pi.clauses[c1].lit.size() > pi.clauses[c2.F].lit.size() && pi.canSubsume(c2.F, c1) && pi.clauses[c2.F].isHard()) {
					if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2.F].lit)) {
						toRemove.push_back(c1);
					}
				}
				else if (pi.clauses[c2.F].lit.size() > pi.clauses[c1].lit.size() && pi.canSubsume(c1, c2.F) && pi.clauses[c1].isHard()) {
					if (isSubsumed(pi.clauses[c2.F].lit, pi.clauses[c1].lit)) {
						toRemove.push_back(c2.F);
					}
				}
			}
		}
	}
}


void Preprocessor::trySEAmsLex(vector<int>& clauses, vector<int>& toRemove) {
	vector<int> cand = amsLex.amsLexSE(clauses);
	for (int c1 : cand) {
		for (int c2 : clauses) {
			if (c1 == c2) continue;
			if (!pi.canSubsume1(c1, c2)) continue;
			if (pi.clauses[c1].lit.size() == pi.clauses[c2].lit.size() && pi.canSubsume(c1, c2) && pi.canSubsume(c2, c1) && c1 > c2) {
				if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2].lit)) {
					if (pi.clauses[c1].isHard() && pi.clauses[c2].isHard()) {
						toRemove.push_back(c2);
					}
					else if(pi.clauses[c1].isHard()) {
						toRemove.push_back(c2);
					}
					else if(pi.clauses[c2].isHard()) {
						toRemove.push_back(c1);
					}
					else {
						pi.clauses[c1].weight += pi.clauses[c2].weight;
						pi.clauses[c2].weight = 0;
						toRemove.push_back(c2);
					}
				}
			}
			else if (pi.clauses[c1].lit.size() > pi.clauses[c2].lit.size() && pi.canSubsume(c2, c1) && pi.clauses[c2].isHard()) {
				if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2].lit)) {
					toRemove.push_back(c1);
				}
			}
			else if (pi.clauses[c2].lit.size() > pi.clauses[c1].lit.size() && pi.canSubsume(c1, c2) && pi.clauses[c1].isHard()) {
				if (isSubsumed(pi.clauses[c2].lit, pi.clauses[c1].lit)) {
					toRemove.push_back(c2);
				}
			}
		}
	}
}

void Preprocessor::trySE(vector<int>& clauses, vector<int>& toRemove) {
	for (int c1 : clauses) {
		for (int c2 : clauses) {
			if (c1 <= c2) continue;
			if (!pi.canSubsume1(c1, c2)) continue;
			if (pi.clauses[c1].lit.size() == pi.clauses[c2].lit.size() && pi.canSubsume(c1, c2) && pi.canSubsume(c2, c1)) {
				if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2].lit)) {
					if (pi.clauses[c1].isHard() && pi.clauses[c2].isHard()) {
						toRemove.push_back(c2);
					}
					else if(pi.clauses[c1].isHard()) {
						toRemove.push_back(c2);
					}
					else if(pi.clauses[c2].isHard()) {
						toRemove.push_back(c1);
					}
					else {
						pi.clauses[c1].weight += pi.clauses[c2].weight;
						pi.clauses[c2].weight = 0;
						toRemove.push_back(c2);
					}
				}
			}
			else if (pi.clauses[c1].lit.size() > pi.clauses[c2].lit.size() && pi.canSubsume(c2, c1) && pi.clauses[c2].isHard()) {
				if (isSubsumed(pi.clauses[c1].lit, pi.clauses[c2].lit)) {
					toRemove.push_back(c1);
				}
			}
			else if (pi.clauses[c2].lit.size() > pi.clauses[c1].lit.size() && pi.canSubsume(c1, c2) && pi.clauses[c1].isHard()) {
				if (isSubsumed(pi.clauses[c2].lit, pi.clauses[c1].lit)) {
					toRemove.push_back(c2);
				}
			}
		}
	}
}

void Preprocessor::trySEgen(int lit, vector<int>& toRemove) {
	if ((int)pi.litClauses[lit].size() <= 64) { // magic constant
		trySE(pi.litClauses[lit], toRemove);
	}
	else {
		int avg = 0;
		for (int c : pi.litClauses[lit]) {
			avg += pi.clauses[c].lit.size();
		}
		avg /= (int)pi.litClauses[lit].size();
		
		if (avg <= 10) { // magic constant
			trySEHash(pi.litClauses[lit], lit, toRemove);
		}
		else if (avg >= 30) { // magic constant
			trySEAmsLex(pi.litClauses[lit], toRemove);
		}
		else {
			if (avg * (int)pi.litClauses[lit].size() > 10000) { // magic constant
				trySEAmsLex(pi.litClauses[lit], toRemove);
			}
			else {
				trySEHash(pi.litClauses[lit], lit, toRemove);
			}
		}
	}
}

int Preprocessor::doSE() {
	rLog.startTechnique(Log::Technique::SE);
	if (!rLog.requestTime(Log::Technique::SE)) {
		rLog.stopTechnique(Log::Technique::SE);
		return 0;
	}
	
	vector<int> checkLit = pi.tl.getModLiterals("SE");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int lit1, int lit2) {
			return pi.litClauses[lit1].size() < pi.litClauses[lit2].size();
		};
		sort(checkLit.begin(), checkLit.end(), cmp);
	}
	bool skip = false;
	vector<int> toRemove;
	if (skipTechnique > 0 && (int)checkLit.size() >= skipTechnique*4) { // magic constant
		for (int tc = 0; tc < skipTechnique; tc++) {
			if (!rLog.requestTime(Log::Technique::SE)) break;
			int lit = checkLit[getRand(0, (int)checkLit.size() - 1)];
			trySEgen(lit, toRemove);
		}
		if (toRemove.size() == 0) {
			skip = true;
			log("SE skipped");
		}
	}
	if (!skip) {
		for (int lit : checkLit) {
			if (!rLog.requestTime(Log::Technique::SE)) break;
			trySEgen(lit, toRemove);
		}
	}
	int removed = 0;
	for (int c : toRemove) {
		if (!pi.isClauseRemoved(c)) {
			pi.removeClause(c);
			removed++;
		}
	}
	rLog.removeClause(removed);
	log(removed, " clauses removed by SE");
	rLog.stopTechnique(Log::Technique::SE);
	return removed;
}

void Preprocessor::doSE2() {
	for (int lit = 0; lit < 2*pi.vars; lit++) {
		if (trySESlow(lit) != 0) {
			print("fail SE ", litToDimacs(lit));
			abort();
		}
	}
}