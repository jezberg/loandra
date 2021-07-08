// Supposes that the clauses are in sorted order
// Can we do SSR such that var is removed from c2?
bool Preprocessor::canSSR(int var, const Clause& c1, const Clause& c2) {
	if (c1.lit.size() > c2.lit.size()) return false;
	int i2 = 0;
	for (int i = 0; i < (int)c1.lit.size(); i++) {
		while (i2 < (int)c2.lit.size() && c2.lit[i2] < c1.lit[i]) {
			i2++;
		}
		if (litVariable(c1.lit[i]) != var && (i2 >= (int)c2.lit.size() || c2.lit[i2] != c1.lit[i])) return false;
	}
	return true;
}

bool Preprocessor::SSRC(int c1, int c2, int var) {
	uint64_t bmP = ((uint64_t)1 << (uint64_t)(posLit(var)&63));
	uint64_t bmN = ((uint64_t)1 << (uint64_t)(negLit(var)&63));
	
	bool canP = false;
	bool canN = false;
	
	if (((pi.clauses[c2].hash ^ bmN) | pi.clauses[c1].hash) == pi.clauses[c1].hash) canP = canSSR(var, pi.clauses[c2], pi.clauses[c1]);
	if (((pi.clauses[c1].hash ^ bmP) | pi.clauses[c2].hash) == pi.clauses[c2].hash) canN = canSSR(var, pi.clauses[c1], pi.clauses[c2]);
	if (pi.isClauseRemoved(c1)) return false;
	if (pi.isClauseRemoved(c2)) return false;
	if (canP && canN) {
		pi.removeLiteralFromClause(posLit(var), c1);
		pi.removeClause(c2);
		rLog.removeLiteral(1);
		rLog.removeClause(1);
		return true;
	}
	else if (canP) {
		pi.removeLiteralFromClause(posLit(var), c1);
		rLog.removeLiteral(1);
		return true;
	}
	else if (canN) {
		pi.removeLiteralFromClause(negLit(var), c2);
		rLog.removeLiteral(1);
		return true;
	}
	return false;
}

int Preprocessor::trySSRAmsLex(int var) {
	pair<vector<int>, vector<int> > cand = amsLex.amsLexSSR(pi.litClauses[posLit(var)], pi.litClauses[negLit(var)], var);
	int removed = 0;
	for (int c1 : cand.F) {
		if (pi.isClauseRemoved(c1)) continue;
		for (int c2 : pi.litClauses[negLit(var)]) {
			if (pi.isClauseRemoved(c2)) continue;
			if (!binary_search(pi.clauses[c1].lit.begin(), pi.clauses[c1].lit.end(), posLit(var))) continue;
			if (!binary_search(pi.clauses[c2].lit.begin(), pi.clauses[c2].lit.end(), negLit(var))) continue;
			if (SSRC(c1, c2, var)) {
				removed++;
			}
		}
	}
	for (int c1 : pi.litClauses[posLit(var)]) {
		if (pi.isClauseRemoved(c1)) continue;
		for (int c2 : cand.S) {
			if (pi.isClauseRemoved(c2)) continue;
			if (!binary_search(pi.clauses[c1].lit.begin(), pi.clauses[c1].lit.end(), posLit(var))) continue;
			if (!binary_search(pi.clauses[c2].lit.begin(), pi.clauses[c2].lit.end(), negLit(var))) continue;
			if (SSRC(c1, c2, var)) {
				removed++;
			}
		}
	}
	return removed;
}

int Preprocessor::trySSRHash(int var) {
	uint64_t bmP = ((uint64_t)1 << (uint64_t)(posLit(var)&63));
	uint64_t bmN = ((uint64_t)1 << (uint64_t)(negLit(var)&63));
	
	vector<int>& pc = pi.litClauses[posLit(var)];
	vector<int>& nc = pi.litClauses[negLit(var)];
	uint64_t k = 1;
	while ((1<<k) < max((int)pc.size(), (int)nc.size())) k++;
	
	vector<vector<pair<int, uint64_t> > > hp(1<<k);
	vector<vector<pair<int, uint64_t> > > hn(1<<k);
	
	for (int c : pc) {
		// do UP to avoid special case
		if (pi.clauses[c].lit.size() == 1){
			return setVariable(var, true);
		}
		uint64_t h = 0;
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) != var) h |= ((uint64_t)1 << (uint64_t)(l%k));
		}
		hp[h].push_back({c, pi.clauses[c].hash});
	}
	
	for (int c : nc) {
		// do UP to avoid special case
		if (pi.clauses[c].lit.size() == 1){
			return setVariable(var, false);
		}
		uint64_t h = 0;
		for (int l : pi.clauses[c].lit) {
			if (litVariable(l) != var) h |= ((uint64_t)1 << (uint64_t)(l%k));
		}
		hn[h].push_back({c, pi.clauses[c].hash});
	}
	
	int removed = 0;
	bool f  = true;
	while (f) {
		f = false;
		for (int c1 : pi.litClauses[posLit(var)]) {
			int64_t h = 0;
			for (int l : pi.clauses[c1].lit) {
				if (litVariable(l) != var) h |= ((int64_t)1 << (int64_t)(l%k));
			}
			for (int64_t sub = 0; (sub = (sub - h) & h);) {
				for (auto c2 : hn[sub]) {
					if ((((pi.clauses[c1].hash ^ bmP) | c2.S) != c2.S) && (((c2.S ^ bmN) | pi.clauses[c1].hash) != pi.clauses[c1].hash)) continue;
					if ((c2.S & bmN) == 0) continue;
					if (pi.isClauseRemoved(c2.F)) continue;
					if (!binary_search(pi.clauses[c2.F].lit.begin(), pi.clauses[c2.F].lit.end(), negLit(var))) continue;
					if (SSRC(c1, c2.F, var)) {
						f = true;
						removed++;
						break;
					}
				}
				if (f) break;
			}
			if (f) break;
		}
	}
	f = true;
	while (f) {
		f = false;
		for (int c2 : pi.litClauses[negLit(var)]) {
			int64_t h = 0;
			for (int l : pi.clauses[c2].lit) {
				if (litVariable(l) != var) h |= ((int64_t)1 << (int64_t)(l%k));
			}
			for (int64_t sub = 0; (sub = (sub - h) & h);) {
				for (auto c1 : hp[sub]) {
					if ((((c1.S ^ bmP) | pi.clauses[c2].hash) != pi.clauses[c2].hash) && (((pi.clauses[c2].hash ^ bmN) | c1.S) != c1.S)) continue;
					if ((c1.S & bmP) == 0) continue;
					if (pi.isClauseRemoved(c1.F)) continue;
					if (!binary_search(pi.clauses[c1.F].lit.begin(), pi.clauses[c1.F].lit.end(), posLit(var))) continue;
					if (SSRC(c1.F, c2, var)) {
						f = true;
						removed++;
						break;
					}
				}
				if (f) break;
			}
			if (f) break;
		}
	}
	return removed;
}

int Preprocessor::trySSR2(int var) {
	int removed = 0;
	
	uint64_t bmP = ((uint64_t)1 << (uint64_t)(posLit(var)&63));
	uint64_t bmN = ((uint64_t)1 << (uint64_t)(negLit(var)&63));
	
	bool f = true;
	while (f) {
		f = false;
		for (int c1 : pi.litClauses[posLit(var)]) {
			for (int c2 : pi.litClauses[negLit(var)]) {
				if ((((pi.clauses[c1].hash ^ bmP) | pi.clauses[c2].hash) != pi.clauses[c2].hash) && (((pi.clauses[c2].hash ^ bmN) | pi.clauses[c1].hash) != pi.clauses[c1].hash)) continue;
				
				if (SSRC(c1, c2, var)) {
					f = true;
					removed++;
					break;
				}
			}
			if (f) break;
		}
	}
	return removed;
}

// Trivial implementation
int Preprocessor::trySSR(int var) {
	int removed = 0;
	bool f = true;
	while (f) {
		f = false;
		for (int i = 0; i < (int)pi.litClauses[posLit(var)].size(); i++) {
			for (int ii = 0; ii < (int)pi.litClauses[negLit(var)].size(); ii++) {
				Clause& cp = pi.clauses[pi.litClauses[posLit(var)][i]];
				Clause& cn = pi.clauses[pi.litClauses[negLit(var)][ii]];
				
				bool canP = canSSR(var, cn, cp);
				bool canN = canSSR(var, cp, cn);
				
				if (canP && canN) {
					pi.removeLiteralFromClause(posLit(var), pi.litClauses[posLit(var)][i]);
					rLog.removeLiteral(1);
					pi.removeClause(pi.litClauses[negLit(var)][ii]);
					rLog.removeClause(1);
					removed++;
					f = true;
					break;
				}
				else if (canP) {
					pi.removeLiteralFromClause(posLit(var), pi.litClauses[posLit(var)][i]);
					rLog.removeLiteral(1);
					removed++;
					f = true;
					break;
				}
				else if (canN) {
					pi.removeLiteralFromClause(negLit(var), pi.litClauses[negLit(var)][ii]);
					rLog.removeLiteral(1);
					removed++;
					f = true;
					break;
				}
			}
			if (f) break;
		}
	}
	return removed;
}

int Preprocessor::trySSRgen(int var) {
	int eliminated = 0;
	if (pi.isLabel[var] == 0 && pi.litClauses[posLit(var)].size() > 0 && pi.litClauses[negLit(var)].size() > 0) {
		uint64_t n = (uint64_t)pi.litClauses[posLit(var)].size()*(uint64_t)pi.litClauses[negLit(var)].size();
		if (n <= 50000) { // magic constant
			eliminated += trySSR2(var);
		}
		else {
			int avg = 0;
			for (int c : pi.litClauses[posLit(var)]) {
				avg += pi.clauses[c].lit.size();
			}
			for (int c : pi.litClauses[negLit(var)]) {
				avg += pi.clauses[c].lit.size();
			}
			avg /= (int)(pi.litClauses[posLit(var)].size() + pi.litClauses[negLit(var)].size());
			
			if (avg <= 4) { // magic constant
				eliminated += trySSRHash(var);
			}
			else if (avg >= 8) { // magic constant
				eliminated += trySSRAmsLex(var);
			}
			else {
				if (avg * (int)(pi.litClauses[posLit(var)].size() + pi.litClauses[negLit(var)].size()) > 10000) { // magic constant
					eliminated += trySSRAmsLex(var);
				}
				else {
					eliminated += trySSRHash(var);
				}
			}
		}
	}
	return eliminated;
}

int Preprocessor::doSSR() {
	rLog.startTechnique(Log::Technique::SSR);
	if (!rLog.requestTime(Log::Technique::SSR)) {
		rLog.stopTechnique(Log::Technique::SSR);
		return 0;
	}
	int eliminated = 0;
	vector<int> checkVar = pi.tl.getModVariables("SSR");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[negLit(var1)].size() + pi.litClauses[posLit(var1)].size() < pi.litClauses[negLit(var2)].size() + pi.litClauses[posLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	bool skip = false;
	if (skipTechnique > 0 && (int)checkVar.size() >= skipTechnique * 4) { // magic constant
		for (int tc = 0; tc < skipTechnique; tc++) {
			if (!rLog.requestTime(Log::Technique::SSR)) break;
			int var = checkVar[getRand(0, (int)checkVar.size() -1)];
			eliminated += trySSRgen(var);
		}
		if (eliminated == 0) {
			skip = true;
			log("SSR skipped");
		}
	}
	if (!skip) {
		for (int var : checkVar) {
			if (!rLog.requestTime(Log::Technique::SSR)) break;
			eliminated += trySSRgen(var);
		}
	}
	log(eliminated, " eliminated by SSR");
	rLog.stopTechnique(Log::Technique::SSR);
	return eliminated;
}

void Preprocessor::doSSR2() {
	for (int var = 0; var < pi.vars; var++) {
		if (pi.isLabel[var] == 0 && pi.litClauses[posLit(var)].size() > 0 && pi.litClauses[negLit(var)].size() > 0) {
			if (trySSR(var) != 0){
				print("fail SSR ", var + 1);
				abort();
			}
		}
	}
}