void Preprocessor::GSLEBT(int i, uint64_t w, vector<int>& sel, vector<uint64_t>& weights, vector<vector<int> >& hs, bool& found, uint64_t& itLim) {
	if (i == (int)hs.size()) {
		found = true;
		return;
	}
	if (itLim == 0) {
		return;
	}
	itLim--;
	bool f = false;
	for (int l : hs[i]) {
		if (sel[l]) {
			f = true;
			break;
		}
	}
	if (f) {
		GSLEBT(i + 1, w, sel, weights, hs, found, itLim);
	}
	else {
		for (int l : hs[i]) {
			if (weights[l] > w) continue;
			sel[l] = 1;
			GSLEBT(i + 1, w - weights[l], sel, weights, hs, found, itLim);
			sel[l] = 0;
		}
	}
}

bool Preprocessor::GSLEtryBackTrack(vector<vector<int> >& hs, vector<uint64_t>& weights, uint64_t w, uint64_t itLim) {
	if (hs.size() == 0) return true;
	bool found = false;
	vector<int> sel(weights.size());
	found = false;
	int tries = 10; // magic constant
	for (int l : hs[0]) {
		if (weights[l] > w) continue;
		uint64_t lm = itLim;
		sel[l] = 1;
		GSLEBT(1, w - weights[l], sel, weights, hs, found, lm);
		sel[l] = 0;
		if (found) break;
		tries--;
		if (tries < 0) break;
	}
	return found;
}

int Preprocessor::tryGSLE(int lb) {
	uint64_t lw = pi.labelWeight(litVariable(lb));
	for (int c : pi.litClauses[lb]) {
		bool f = false;
		for (int l : pi.clauses[c].lit) {
			if (l != lb && pi.isLabel[litVariable(l)] && pi.labelWeight(litVariable(l)) <= lw) {
				f = true;
				break;
			}
		}
		if (!f) {
			return 0;
		}
	}
	vector<vector<int> > hs;
	set<int> sel;
	uint64_t selW = 0;
	for (int c : pi.litClauses[lb]) {
		vector<int> s;
		for (int l : pi.clauses[c].lit) {
			if (l != lb && pi.isLabel[litVariable(l)] && pi.labelWeight(litVariable(l)) <= lw) {
				s.push_back(litVariable(l));
			}
		}
		hs.push_back(s);
	}
	for (auto& s : hs) {
		if (s.size() == 1) {
			if (sel.count(s[0]) == 0) {
				sel.insert(s[0]);
				selW += pi.labelWeight(s[0]);
			}
		}
	}
	if (selW > lw) return 0;
	vector<vector<int> > hss;
	uint64_t hsl = 0;
	for (auto& s : hs) {
		if (s.size() > 1) {
			bool f = false;
			for (int l : s) {
				if (sel.count(l) > 0) {
					f = true;
					break;
				}
			}
			if (!f) {
				hss.push_back(s);
				hsl += (uint64_t)(int)s.size();
			}
		}
	}
	auto cmp = [](const vector<int>& a, const vector<int>& b) {
		return a.size() < b.size();
	};
	sort(hss.begin(), hss.end(), cmp);
	// Do coordinate compression
	vector<int> cc;
	for (auto& se : hss) {
		for (int l : se) {
			cc.push_back(l);
		}
	}
	sort(cc.begin(), cc.end());
	cc.erase(unique(cc.begin(), cc.end()), cc.end());
	vector<uint64_t> laWeights(cc.size());
	for (int i = 0; i < (int)cc.size(); i++) {
		laWeights[i] = pi.labelWeight(cc[i]);
	}
	for (auto& se : hss) {
		for (int& l : se) {
			l = lower_bound(cc.begin(), cc.end(), l) - cc.begin();
		}
	}
	vector<int> cnt(cc.size());
	for (int i = (int)hss.size() - 1; i >= 0; i--) {
		for (int l : hss[i]) {
			cnt[l]++;
		}
		auto cmp2 = [&](int a, int b) {// the log-approximation heuristic
			return (double)laWeights[a]/(double)cnt[a] < (double)laWeights[b]/(double)cnt[b];
		};
		sort(hss[i].begin(), hss[i].end(), cmp2);
	}
	bool ok = false;
	hsl = min(hsl, (uint64_t)hss.size()*10); // magic constant
	if (GSLEtryBackTrack(hss, laWeights, lw - selW, hsl)) {
		ok = true;
	}
	if (ok) {
		int rmClauses = 0;
		if (pi.isLabel[litVariable(lb)] == VAR_TRUE) {
			rmClauses = setVariable(litVariable(lb), true);
		}
		else {
			rmClauses = setVariable(litVariable(lb), false);
		}
		assert(pi.isVarRemoved(litVariable(lb)));
		rLog.removeClause(rmClauses);
		rLog.removeLabel(1);
		return 1;
	}
	return 0;
}

int Preprocessor::doGSLE() {
	rLog.startTechnique(Log::Technique::GSLE);
	int removed = 0;
	if (!rLog.requestTime(Log::Technique::GSLE)) {
		rLog.stopTechnique(Log::Technique::GSLE);
		return 0;
	}
	vector<int> checkVar = pi.tl.getTouchedVariables("GSLE");
	if (rLog.isTimeLimit()) {
		auto cmp = [&](int var1, int var2) {
			return pi.litClauses[negLit(var1)].size() + pi.litClauses[posLit(var1)].size() < pi.litClauses[negLit(var2)].size() + pi.litClauses[posLit(var2)].size();
		};
		sort(checkVar.begin(), checkVar.end(), cmp);
	}
	bool skip = false;
	if (skipTechnique > 0 && (int)checkVar.size() >= skipTechnique * 4) { // magic constant
		for (int tc = 0; tc < skipTechnique; tc++) {
			if (!rLog.requestTime(Log::Technique::GSLE)) break;
			int var = checkVar[getRand(0, (int)checkVar.size() - 1)];
			if (pi.isLabel[var] == VAR_UNDEFINED) continue;
			if (pi.isVarRemoved(var)) continue;
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
			
			if (pi.isLabel[var] == VAR_TRUE) {
				removed += tryGSLE(negLit(var));
			}
			else {
				removed += tryGSLE(posLit(var));
			}
		}
		if (removed == 0) {
			skip = true;
			log("GSLE skipped");
		}
	}
	if (!skip) {
		for (int var : checkVar) {
			if (pi.isLabel[var] == VAR_UNDEFINED) continue;
			if (pi.isVarRemoved(var)) continue;
			if (!rLog.requestTime(Log::Technique::GSLE)) break;
			
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
			
			if (pi.isLabel[var] == VAR_TRUE) {
				removed += tryGSLE(negLit(var));
			}
			else {
				removed += tryGSLE(posLit(var));
			}
		}
	}
	
	log(removed, " labels removed by GSLE");
	rLog.stopTechnique(Log::Technique::GSLE);
	return removed;
}

void Preprocessor::doGSLE2() {
	vector<int> lbs;
	for (int var = 0; var < pi.vars; var++) {
		if (pi.isLabel[var] != VAR_UNDEFINED && !pi.isVarRemoved(var)) {
			lbs.push_back(var);
		}
	}
	for (int var : lbs) {
		if (pi.isLabel[var] == VAR_TRUE) {
			if (tryGSLE(negLit(var))) {
				print("fail GSLE");
				print(var + 1);
				abort();
			}
		}
		else {
			if (tryGSLE(posLit(var))) {
				print("fail GSLE");
				print(var + 1);
				abort();
			}
		}
	}
}