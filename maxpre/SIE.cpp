bool Preprocessor::SIErndCheck(int litX, int litY) {
	if (pi.litClauses[litNegation(litX)].size() == 0) return true;
	int c = pi.litClauses[litNegation(litX)].back();
	if (pi.clauses[c].lit.size() > 2) { // magic constant
		if (!binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litNegation(litY))) {
			return false;
		}
	}
	else {
		bool f = false;
		for (int l : pi.clauses[c].lit) {
			if (l == litNegation(litY)) {
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

int Preprocessor::try2SIE(int litX, int litY) {
	for (int c : pi.litClauses[litNegation(litX)]) {
		if (pi.clauses[c].lit.size() > 10) { // magic constant
			if (!binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litNegation(litY))) {
				return 0;
			}
		}
		else {
			bool f = false;
			for (int l : pi.clauses[c].lit) {
				if (l == litNegation(litY)) {
					f = true;
					break;
				}
			}
			if (!f) {
				return 0;
			}
		}
	}
	for (int c : pi.litClauses[litY]) {
		if (!binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litX)) {
			return 0;
		}
	}
	int rmC = 0;
	rmC += setVariable(litVariable(litX), litValue(litX));
	rmC += setVariable(litVariable(litY), !litValue(litY));
	rLog.removeClause(rmC);
	rLog.removeVariable(2);
	assert(pi.isVarRemoved(litVariable(litX)));
	assert(pi.isVarRemoved(litVariable(litY)));
	return rmC;
}

// lit is l_y
int Preprocessor::trySIE(int lit) {
	if (pi.isVarRemoved(litVariable(lit))) return 0;
	if (pi.litClauses[lit].size() == 0) {
		int rm = setVariable(litVariable(lit), !litValue(lit));
		return rm;
	}
	if (pi.litClauses[lit].size() == 1){
		int c = pi.litClauses[lit][0];
		for (int l : pi.clauses[c].lit) {
			if (l != lit) {
				int rm = try2SIE(l, lit);
				if (rm > 0) {
					return rm;
				}
			}
		}
		return 0;
	}
	unsigned mi1 = 1e9;
	int mi1c = -1;
	unsigned mi2 = 1e9;
	int mi2c = -1;
	for (int c : pi.litClauses[lit]) {
		if (pi.clauses[c].lit.size() < mi1) {
			mi2 = mi1;
			mi2c = mi1c;
			mi1 = pi.clauses[c].lit.size();
			mi1c = c;
		}
		else if(pi.clauses[c].lit.size() < mi2) {
			mi2 = pi.clauses[c].lit.size();
			mi2c = c;
		}
	}
	assert(mi1c != mi2c);
	assert(mi1c >= 0 && mi2c >= 0);
	unsigned i1 = 0;
	unsigned i2 = 0;
	vector<int> cands;
	while (i1 < pi.clauses[mi1c].lit.size() && i2 < pi.clauses[mi2c].lit.size()) {
		if (pi.clauses[mi1c].lit[i1] == pi.clauses[mi2c].lit[i2]) {
			if (pi.clauses[mi1c].lit[i1] != lit) {
				cands.push_back(pi.clauses[mi1c].lit[i1]);
			}
			i1++;
			i2++;
		}
		else if(pi.clauses[mi1c].lit[i1] < pi.clauses[mi2c].lit[i2]) {
			i1++;
		}
		else {
			i2++;
		}
	}
	vector<int> newCands;
	for (int l : cands) {
		if (SIErndCheck(l, lit)) newCands.push_back(l);
	}
	cands = newCands;
	for (int c : pi.litClauses[lit]) {
		vector<int> nCands;
		nCands.reserve(cands.size());
		int ci = 0;
		for (int l : pi.clauses[c].lit) {
			while (ci < (int)cands.size() && cands[ci] < l) ci++;
			if (ci == (int)cands.size()) break;
			if (l == cands[ci]) nCands.push_back(l);
		}
		cands = nCands;
		if (cands.size() < 2) break;
	}
	for (int l : cands) {
		int rm = try2SIE(l, lit);
		if (rm > 0) {
			return rm;
		}
	}
	return 0;
}

int Preprocessor::doSIE() {
	rLog.startTechnique(Log::Technique::SIE);
	int removed = 0;
	vector<int> checkLit = pi.tl.getTouchedLiterals("SIE");
	for (int lit : checkLit) {
		removed += trySIE(lit);
	}
	log(removed, " clauses removed by SIE");
	rLog.stopTechnique(Log::Technique::SIE);
	return removed;
}

void Preprocessor::doSIE2() {
	for (int litX = 0; litX < pi.vars*2; litX++) {
		for (int litY = 0; litY < pi.vars*2; litY++) {
			if (litVariable(litX) != litVariable(litY) && pi.isLabel[litVariable(litX)] == VAR_UNDEFINED && pi.isLabel[litVariable(litY)] == VAR_UNDEFINED) {
				if (!pi.isVarRemoved(litVariable(litX)) && !pi.isVarRemoved(litVariable(litY))) {
					if (try2SIE(litX, litY)) {
						print("fail SIE ", litToDimacs(litX), " ", litToDimacs(litY));
						print("fail SIE ");
						abort();
					}
				}
			}
		}
	}
}