
int Preprocessor::tryFLP(vector<int> fLit, int clause) {
	bool is = false;
	for (int lit : fLit) {
		for (int c : pi.litClauses[lit]) {
			if (pi.clauses[c].lit.size() < 3) {
				is = true;
				break;
			}
		}
	}
	if (!is) return 0;
	queue<int> failQ;
	for (int lit : fLit) {
		failQ.push(lit);
	}
	vector<pair<int, int> > litRemoved;
	bool conflict = false;
	unordered_set<int> satClauses;
	unordered_set<int> used;
	vector<int> binCores;
	while (!failQ.empty()) {
		int lit = failQ.front();
		failQ.pop();
		if (used.count(lit)) continue;
		used.insert(lit);
		for (int c : pi.litClauses[litNegation(lit)]) {
			satClauses.insert(c);
		}
		vector<int> cs = pi.litClauses[lit];
		for (int c : cs) {
			if (satClauses.count(c)) continue;
			if (!pi.clauses[c].isHard()) {
				if (fLit.size() == 1) binCores.push_back(litNegation(pi.clauses[c].lit[0]));
				continue;
			}
			if (!pi.clauses[c].isHard()) continue;
			if (pi.clauses[c].lit.size() == 1) {
				conflict = true;
				break;
			}
			pi.removeLiteralFromClause(lit, c, false);
			litRemoved.push_back({lit, c});
			if (pi.clauses[c].lit.size() == 1) {
				failQ.push(litNegation(pi.clauses[c].lit[0]));
			}
		}
		if (conflict) break;
	}
	for (auto t : litRemoved) {
		pi.addLiteralToClause(t.F, t.S, false);
	}
	if (conflict) {
		if (fLit.size() == 1) {
			int rmClauses;
			if (litValue(fLit[0]) == true) {
				rmClauses = setVariable(litVariable(fLit[0]), true);
			}
			else {
				rmClauses = setVariable(litVariable(fLit[0]), false);
			}
			assert(pi.isVarRemoved(litVariable(fLit[0])));
			rLog.removeClause(rmClauses);
			if (pi.isLabel[litVariable(fLit[0])]) rLog.removeLabel(1);
			else rLog.removeVariable(1);
		}
		else {
			vector<int> teClause = pi.clauses[clause].lit;
			for (int lit : teClause) {
				if (!pi.isLabel[litVariable(lit)]) {
					pi.removeLiteralFromClause(lit, clause);
					rLog.removeLiteral(1);
				}
			}
		}
		return 1;
	}
	return 0;
}

int Preprocessor::doFLP() {
	rLog.startTechnique(Log::Technique::FLP);
	if (!rLog.requestTime(Log::Technique::FLP)) {
		rLog.stopTechnique(Log::Technique::FLP);
		return 0;
	}
	int removed = 0;
	for (int c = 0; c < (int)pi.clauses.size(); c++) {
		if (!rLog.requestTime(Log::Technique::FLP)) break;
		if (pi.isClauseRemoved(c)) continue;
		if (!pi.clauses[c].isHard()) continue;
		vector<int> lbs;
		bool f = false;
		for (int lit : pi.clauses[c].lit) {
			if (pi.isLabel[litVariable(lit)]) lbs.push_back(lit);
			else f = true;
		}
		if (lbs.size() > 0 && f) removed += tryFLP(lbs, c);
	}
	log(removed, " cores found by FLP");
	rLog.stopTechnique(Log::Technique::FLP);
	return removed;
}