void Preprocessor::tryLSBCE(int lit, unordered_set<int>& deletedClauses, unordered_set<int>& touchedList, vector<pair<int, int> >& blockedClauses) {
	for (int c1i : pi.litClauses[lit]) {
		if (deletedClauses.count(c1i)) continue;
		const vector<int>& c1 = pi.clauses[c1i].lit;
		bool f = false;
		for (int c2i : pi.litClauses[litNegation(lit)]) {
			if (deletedClauses.count(c2i)) continue;
			const vector<int>& c2 = pi.clauses[c2i].lit;
			bool ff = false;
			for (int c1l : c1) {
				for (int c2l : c2) {
					if (litNegation(c1l) == c2l && c1l != lit) {
						ff = true;
						break;
					}
				}
				if (ff) break;
			}
			if (!ff) {
				f = true;
				break;
			}
		}
		if (!f) {
			deletedClauses.insert(c1i);
			blockedClauses.push_back({c1i, lit});
			for (int l : c1) {
				touchedList.insert(litVariable(l));
			}
		}
	}
}

int Preprocessor::tryLS(int lbl) {
	assert(pi.isLabel[lbl]);
	vector<int> tc;
	if (pi.isLabel[lbl] == VAR_TRUE) {
		tc = pi.litClauses[negLit(lbl)];
	}
	else {
		tc = pi.litClauses[posLit(lbl)];
	}
	if (tc.size() == 0) return 0;
	vector<pair<int, int> > blockedClauses;
	unordered_set<int> touchedList;
	unordered_set<int> deletedClauses;
	for (int t : tc) {
		deletedClauses.insert(t);
		for(int l : pi.clauses[t].lit) {
			touchedList.insert(litVariable(l));
		}
	}
	while (touchedList.size() > 0) {
		int var = *touchedList.begin();
		touchedList.erase(var);
		if (pi.isLabel[var]) continue;
		tryLSBCE(posLit(var), deletedClauses, touchedList, blockedClauses);
		tryLSBCE(negLit(var), deletedClauses, touchedList, blockedClauses);
	}
	for (auto c : blockedClauses) {
		if (pi.isLabel[lbl] == VAR_TRUE) trace.LS(negLit(lbl), c.S, pi.clauses[c.F].lit);
		else trace.LS(posLit(lbl), c.S, pi.clauses[c.F].lit);
		
		if(pi.isLabel[lbl] == VAR_TRUE) pi.addLiteralToClause(negLit(lbl), c.F);
		else pi.addLiteralToClause(posLit(lbl), c.F);
	}
	rLog.removeLiteral(-(int)blockedClauses.size());
	return blockedClauses.size();
}

int Preprocessor::doLS() {
	while (doBCE());
	rLog.startTechnique(Log::Technique::LS);
	if (!rLog.requestTime(Log::Technique::LS)) {
		rLog.stopTechnique(Log::Technique::LS);
		return 0;
	}
	int added = 0;
	vector<int> checkVar = pi.tl.getTouchedVariables("LS");
	for (int var : checkVar) {
		if (pi.isLabel[var] == VAR_UNDEFINED) continue;
		if (pi.isVarRemoved(var)) continue;
		if (!rLog.requestTime(Log::Technique::LS)) break;
		added += tryLS(var);
	}
	log(added, " clauses labeled by LS");
	rLog.stopTechnique(Log::Technique::LS);
	return added;
}