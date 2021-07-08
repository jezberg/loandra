bool Preprocessor::vSubsumed(vector<int>& v1, vector<int>& v2) {
	for (int a : v1) {
		bool f = false;
		for (int b : v2) {
			if (a == b) {
				f = true;
				break;
			}
		}
		if (!f) {
			return false;
		}
	}
	return true;
}

// Even though this is not optimal implementation this is fast enough
int Preprocessor::trySLESlow(int lb1, int lb2) {
	vector<int> cs1;
	vector<int> cs2;
	uint64_t w1, w2;
	if (pi.isLabel[lb1] == VAR_TRUE) {
		cs1 = pi.litClauses[negLit(lb1)];
		w1 = pi.clauses[pi.litClauses[posLit(lb1)][0]].weight;
	}
	else {
		cs1 = pi.litClauses[posLit(lb1)];
		w1 = pi.clauses[pi.litClauses[negLit(lb1)][0]].weight;
	}
	if (pi.isLabel[lb2] == VAR_FALSE) {
		cs2 = pi.litClauses[posLit(lb2)];
		w2 = pi.clauses[pi.litClauses[negLit(lb2)][0]].weight;
	}
	else {
		cs2 = pi.litClauses[negLit(lb2)];
		w2 = pi.clauses[pi.litClauses[posLit(lb2)][0]].weight;
	}
	assert(w1 < HARDWEIGHT);
	assert(w2 < HARDWEIGHT);
	
	bool s1 = vSubsumed(cs1, cs2);
	bool s2 = vSubsumed(cs2, cs1);
	
	int rmClauses = 0;
	
	if (s1 && w2 <= w1) {
		if (pi.isLabel[lb1] == VAR_TRUE) {
			rmClauses = setVariable(lb1, true);
		}
		else {
			rmClauses = setVariable(lb1, false);
		}
		assert(pi.isVarRemoved(lb1));
		rLog.removeClause(rmClauses);
		rLog.removeLabel(1);
		return 1;
	}
	else if(s2 && w1 <= w2) {
		if (pi.isLabel[lb2] == VAR_TRUE) {
			rmClauses = setVariable(lb2, true);
		}
		else {
			rmClauses = setVariable(lb2, false);
		}
		assert(pi.isVarRemoved(lb2));
		rLog.removeClause(rmClauses);
		rLog.removeLabel(1);
		return 1;
	}
	return 0;
}

int Preprocessor::doSLE() {
	rLog.startTechnique(Log::Technique::SLE);
	int removed = 0;
	if (!rLog.requestTime(Log::Technique::SLE)) {
		rLog.stopTechnique(Log::Technique::SLE);
		return 0;
	}
	vector<int> checkVar = pi.tl.getTouchedVariables("SLE");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[negLit(var1)].size() + pi.litClauses[posLit(var1)].size() < pi.litClauses[negLit(var2)].size() + pi.litClauses[posLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	for (int var : checkVar) {
		if (pi.isLabel[var] == VAR_UNDEFINED) continue;
		if (pi.isVarRemoved(var)) continue;
		if (!rLog.requestTime(Log::Technique::SLE)) break;
		
		if (pi.isLabel[var] == VAR_TRUE && pi.litClauses[negLit(var)].size() == 0){
			setVariable(var, true);
			removed++;
			continue;
		}
		if (pi.isLabel[var] == VAR_FALSE && pi.litClauses[posLit(var)].size() == 0){
			setVariable(var, false);
			removed++;
			continue;
		}
		
		bool f = true;
		while (f) {
			f = false;
			vector<int> clauses;
			if (pi.isLabel[var] == VAR_TRUE) clauses = pi.litClauses[negLit(var)];
			else clauses = pi.litClauses[posLit(var)];
			for (int c : clauses) {
				for (int l : pi.clauses[c].lit) {
					if (pi.isLabel[litVariable(l)] && litVariable(l) != var) {
						if (trySLESlow(var, litVariable(l))) {
							removed++;
							f = true;
							break;
						}
					}
				}
				if (f) break;
			}
		}
	}
	
	log(removed, " labels removed by SLE");
	rLog.stopTechnique(Log::Technique::SLE);
	return removed;
}

void Preprocessor::doSLE2() {
	vector<int> lbs;
	for (int var = 0; var < pi.vars; var++) {
		if (pi.isLabel[var] != VAR_UNDEFINED && !pi.isVarRemoved(var)) {
			lbs.push_back(var);
		}
	}
	for (int lb1 : lbs) {
		for (int lb2 : lbs) {
			if (lb1 != lb2 && !pi.isVarRemoved(lb1) && !pi.isVarRemoved(lb2)) {
				if (trySLESlow(lb1, lb2)) {
					print("fail SLE");
					print(lb1 + 1, " ", lb2 + 1);
					abort();
				}
			}
		}
	}
}