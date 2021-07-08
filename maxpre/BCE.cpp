int Preprocessor::tryBCE(int lit) {
	vector<int> toRemove;
	for (int c1i : pi.litClauses[lit]) {
		const vector<int>& c1 = pi.clauses[c1i].lit;
		bool f = false;
		for (int c2i : pi.litClauses[litNegation(lit)]) {
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
			toRemove.push_back(c1i);
		}
	}
	for (int c : toRemove) {
		trace.BCE(lit, pi.clauses[c].lit);
		pi.removeClause(c);
	}
	rLog.removeClause((int)toRemove.size());
	return toRemove.size();
}

int Preprocessor::doBCE() {
	rLog.startTechnique(Log::Technique::BCE);
	if (!rLog.requestTime(Log::Technique::BCE)) {
		rLog.stopTechnique(Log::Technique::BCE);
		return 0;
	}
	int removed = 0;
	vector<int> checkVar = pi.tl.getTouchedVariables("BCE");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[posLit(var1)].size() + pi.litClauses[negLit(var1)].size() < pi.litClauses[posLit(var2)].size() + pi.litClauses[negLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	bool skip = false;
	if (skipTechnique > 0 && (int)checkVar.size() >= skipTechnique * 4) {
		for (int tc = 0; tc < skipTechnique; tc++) {
			if (!rLog.requestTime(Log::Technique::BCE)) break;
			int var = checkVar[getRand(0, (int)checkVar.size() - 1)];
			if (pi.isLabel[var] == 0) {
				removed += tryBCE(negLit(var));
				removed += tryBCE(posLit(var));
			}
		}
		if (removed == 0) {
			log("BCE skipped");
			skip = true;
		}
	}
	if (!skip) {
		for (int var : checkVar) {
			if (!rLog.requestTime(Log::Technique::BCE)) break;
			if (pi.isLabel[var] == 0) {
				removed += tryBCE(negLit(var));
				removed += tryBCE(posLit(var));
			}
		}
	}
	log(removed, " removed by BCE");
	rLog.stopTechnique(Log::Technique::BCE);
	return removed;
}

void Preprocessor::doBCE2() {
	for (int lit = 0; lit < 2*pi.vars; lit++) {
		if (pi.isLabel[litVariable(lit)] == 0 && !pi.isVarRemoved(litVariable(lit))) {
			assert(tryBCE(lit) == 0);
		}
	}
}