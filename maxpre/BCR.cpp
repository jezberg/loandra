bool Preprocessor::tryBCR(int c, int l11) {
	int l1 = pi.clauses[c].lit[0];
	int l2 = pi.clauses[c].lit[1];
	if (l2 == l11) swap(l1, l2);
	assert(l1 == l11);
	
	if (pi.litClauses[l1].size() < 2) return false;
	if (pi.litClauses[l2].size() < 2) return false;
	assert(pi.litClauses[litNegation(l1)].size() == 1);
	assert(pi.litClauses[litNegation(l2)].size() == 1);
	
	vector<vector<int> > newClauses;
	
	int sizeLimit = pi.litClauses[l1].size() + pi.litClauses[l2].size();
	
	for (int c1 : pi.litClauses[l1]) {
		if (c1 != c) {
			for (int l : pi.clauses[c1].lit) {
				if (l == l2) return false;
			}
		}
		for (int c2 : pi.litClauses[l2]) {
			if (c1 == c || c2 == c) continue;
			vector<int> nc;
			for (int l : pi.clauses[c1].lit) {
				if (l != l1) nc.push_back(l);
			}
			for (int l : pi.clauses[c2].lit) {
				if (l != l1) nc.push_back(l);
			}
			sort(nc.begin(), nc.end());
			nc.erase(unique(nc.begin(), nc.end()), nc.end());
			bool taut = false;
			for (unsigned i = 1; i < nc.size(); i++) {
				if (litNegation(nc[i]) == nc[i - 1]) {
					taut = true;
					break;
				}
			}
			if (taut) continue;
			if ((int)newClauses.size() >= sizeLimit) return false;
			newClauses.push_back(nc);
		}
	}
	vector<int> toRemove;
	vector<vector<int> > nClauses;
	pi.removeClause(pi.litClauses[litNegation(l1)][0]);
	for (int cc : pi.litClauses[l1]) {
		if (!pi.isClauseRemoved(cc)) {
			pi.removeClause(cc);
			rLog.removeClause(1);
		}
		if (cc != c) {
			nClauses.push_back(vector<int>());
			for (int l : pi.clauses[cc].lit) {
				nClauses.back().push_back(l);
			}
		}
	}
	for (int cc : pi.litClauses[l2]) {
		if (!pi.isClauseRemoved(cc)) {
			pi.removeClause(cc);
			rLog.removeClause(1);
		}
	}
	for (auto& nc : newClauses) {
		pi.addClause(nc);
		rLog.removeClause(-1);
	}
	rLog.removeLabel(1);
	trace.removeWeight(pi.clauses[pi.litClauses[litNegation(l2)][0]].weight);
	trace.BCR(l1, l2, nClauses);
	return true;
}

int bcrc = 0;

int Preprocessor::doBCR() {
	rLog.startTechnique(Log::Technique::BCR);
	int removed = 0;
	if (!rLog.requestTime(Log::Technique::BCR)) {
		rLog.stopTechnique(Log::Technique::BCR);
		return 0;
	}
	vector<int> checkVar = pi.tl.getTouchedVariables("BCR");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[negLit(var1)].size() + pi.litClauses[posLit(var1)].size() < pi.litClauses[negLit(var2)].size() + pi.litClauses[posLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	for (int var : checkVar) {
		if (pi.isLabel[var] == VAR_UNDEFINED) continue;
		if (pi.isVarRemoved(var)) continue;
		if (!rLog.requestTime(Log::Technique::BCR)) break;
		if (pi.isLabel[var] == VAR_TRUE) {
			assert(pi.litClauses[posLit(var)].size() == 1);
			for (int c : pi.litClauses[negLit(var)]) {
				if (pi.clauses[c].lit.size() == 2) {
					if (pi.isLabel[litVariable(pi.clauses[c].lit[0])] && pi.isLabel[litVariable(pi.clauses[c].lit[1])]) {
						if (pi.labelWeight(litVariable(pi.clauses[c].lit[0])) == pi.labelWeight(litVariable(pi.clauses[c].lit[1]))) {
							if (tryBCR(c, negLit(var))) {
								removed++;
								break;
							}
						}
					}
				}
			}
		}
		else if(pi.isLabel[var] == VAR_FALSE) {
			assert(pi.litClauses[negLit(var)].size() == 1);
			for (int c : pi.litClauses[posLit(var)]) {
				if (pi.clauses[c].lit.size() == 2) {
					if (pi.isLabel[litVariable(pi.clauses[c].lit[0])] && pi.isLabel[litVariable(pi.clauses[c].lit[1])]) {
						if (pi.labelWeight(litVariable(pi.clauses[c].lit[0])) == pi.labelWeight(litVariable(pi.clauses[c].lit[1]))) {
							if (tryBCR(c, posLit(var))) {
								removed++;
								break;
							}
						}
					}
				}
			}
		}
		else{
			assert(0);
		}
	}
	rLog.stopTechnique(Log::Technique::BCR);
	log(removed, " labels removed by BCR");
	return removed;
}

void Preprocessor::doBCR2() {
	for (int i = 0; i < pi.vars; i++) {
		if (pi.isLabel[i] && !pi.isVarRemoved(i)) {
			if (pi.isLabel[i] == VAR_TRUE) {
				assert(pi.litClauses[posLit(i)].size() == 1);
				for (int c : pi.litClauses[negLit(i)]) {
					if (pi.clauses[c].lit.size() == 2) {
						if (pi.isLabel[litVariable(pi.clauses[c].lit[0])] && pi.isLabel[litVariable(pi.clauses[c].lit[1])]) {
							if (pi.clauses[pi.litClauses[litNegation(pi.clauses[c].lit[0])][0]].weight == pi.clauses[pi.litClauses[litNegation(pi.clauses[c].lit[1])][0]].weight) {
								if (tryBCR(c, negLit(i))) {
									print("FAIL BCR ", i);
									abort();
								}
							}
						}
					}
				}
			}
			else if(pi.isLabel[i] == VAR_FALSE) {
				assert(pi.litClauses[negLit(i)].size() == 1);
				for (int c : pi.litClauses[posLit(i)]) {
					if (pi.clauses[c].lit.size() == 2) {
						if (pi.isLabel[litVariable(pi.clauses[c].lit[0])] && pi.isLabel[litVariable(pi.clauses[c].lit[1])]) {
							if (pi.clauses[pi.litClauses[litNegation(pi.clauses[c].lit[0])][0]].weight == pi.clauses[pi.litClauses[litNegation(pi.clauses[c].lit[1])][0]].weight) {
								if (tryBCR(c, posLit(i))) {
									print("FAIL BCR ", i);
									abort();
								}
							}
						}
					}
				}
			}
			else{
				assert(0);
			}
		}
	}
}