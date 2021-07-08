// Search all nodes in undirected connected component
void Preprocessor::BIGdfs1(int x, vector<int>& ns) {
	// we do actually bfs to avoid stack overflow
	queue<int> bfs;
	bfs.push(x);
	while (!bfs.empty()) {
		int v = bfs.front();
		bfs.pop();
		if (BIGu[v] == BIGIt) continue;
		BIGu[v] = BIGIt;
		BIGu2[v] = 0;
		BIGu2[litVariable(v)] = 0;
		ns.push_back(v);
		// backward edges
		for (int c : pi.litClauses[v]) {
			if (pi.clauses[c].lit.size() == 2) {
				for (int nx : pi.clauses[c].lit) {
					if (nx != v && pi.isLabel[litVariable(nx)] == VAR_UNDEFINED) {
						bfs.push(litNegation(nx));
					}
				}
			}
		}
		// forward edges
		for (int c : pi.litClauses[litNegation(v)]) {
			if (pi.clauses[c].lit.size() == 2) {
				for (int nx : pi.clauses[c].lit) {
					if (nx != litNegation(v) && pi.isLabel[litVariable(nx)] == VAR_UNDEFINED) {
						bfs.push(nx);
					}
				}
			}
		}
	}
}

// The 2 dfs of Kosaraju
void Preprocessor::BIGdfs2(vector<vector<int> >& g, int x, vector<int>& ns) {
	stack<pair<int, int>, vector<pair<int, int> >  > dfs;
	dfs.push({x, 0});
	while (!dfs.empty()) {
		int v = dfs.top().F;
		int ci = dfs.top().S;
		dfs.pop();
		if (ci == 0 && BIGu2[v] == 1) continue;
		BIGu2[v] = 1;
		if (ci >= (int)g[v].size()) {
			ns.push_back(v);
		}
		else {
			dfs.push({v, ci + 1});
			dfs.push({g[v][ci], 0});
		}
	}
}

void Preprocessor::BIGdfs3(vector<vector<int> >& rg, int x, vector<int>& scc) {
	stack<pair<int, int>, vector<pair<int, int> > > dfs;
	dfs.push({x, 0});
	while (!dfs.empty()) {
		int v = dfs.top().F;
		int ci = dfs.top().S;
		dfs.pop();
		if (ci == 0 && BIGu2[v] == 2) continue;
		BIGu2[v] = 2;
		if (ci == 0) scc.push_back(v);
		if (ci < (int)rg[v].size()) {
			dfs.push({v, ci + 1});
			dfs.push({rg[v][ci], 0});
		}
	}
}

bool Preprocessor::BIGisPath(int from, int to, vector<int>& leIndex, vector<int>& riIndex, vector<int>& up1, vector<int>& up2) {
	int a = up1[from];
	int b = up2[to];
	return (leIndex[a] <= leIndex[to] && riIndex[to] <= riIndex[a]) && (leIndex[b] <= leIndex[from] && riIndex[from] <= riIndex[b]);
}

void Preprocessor::genIndex(vector<vector<int> >& g, vector<vector<int> >& rg, int x, int u1, int u2, int& stamp, vector<int>& le, vector<int>& ri, vector<int>& up1, vector<int>& up2, int order) {
	if (BIGu2[x] == 4) return;
	BIGu2[x] = 4;
	up1[x] = u1;
	up2[x] = u2;
	le[x] = stamp++;
	vector<pair<int, int> > ne;
	ne.reserve(g[x].size() + rg[x].size());
	if (order != 1) {
		for (int nx : g[x]) {
			ne.push_back({nx, 1});
		}
		for (int nx : rg[x]) {
			ne.push_back({nx, 2});
		}
	}
	else {
		for (int nx : rg[x]) {
			ne.push_back({nx, 2});
		}
		for (int nx : g[x]) {
			ne.push_back({nx, 1});
		}
	}
	if (order == 0) {
		random_shuffle(ne.begin(), ne.begin() + g[x].size());
		random_shuffle(ne.begin() + g[x].size(), ne.end());
	}
	else if (order == 1) {
		random_shuffle(ne.begin(), ne.begin() + rg[x].size());
		random_shuffle(ne.begin() + rg[x].size(), ne.end());
	}
	else if (order == 2) {
		random_shuffle(ne.begin(), ne.end());
	}
	else {
		assert(0);
	}
	for (auto nx : ne) {
		if (nx.S == 1) genIndex(g, rg, nx.F, nx.F, u2, stamp, le, ri, up1, up2, order);
		else genIndex(g, rg, nx.F, u1, nx.F, stamp, le, ri, up1, up2, order);
	}
	ri[x] = stamp++;
}

int Preprocessor::tryBIG(int lit, bool doTC) {
	vector<int> cc;
	BIGdfs1(lit, cc);
	for (unsigned i = 0; i < cc.size(); i++) {
		BIGu[cc[i]] = BIGIt - 1;
		BIGid[cc[i]] = (int)i;
		BIGu2[i] = 0;
	}
	vector<vector<int> > g(cc.size());
	vector<vector<int> > rg(cc.size());
	for (int x : cc) {
		// forward edges
		for (int c : pi.litClauses[litNegation(x)]) {
			if (pi.clauses[c].lit.size() == 2) {
				for (int nx : pi.clauses[c].lit) {
					if (nx != litNegation(x) && pi.isLabel[litVariable(nx)] == VAR_UNDEFINED && BIGu[nx] != BIGIt) {
						g[BIGid[x]].push_back(BIGid[nx]);
					}
				}
			}
		}
		// backward edges
		for (int c : pi.litClauses[x]) {
			if (pi.clauses[c].lit.size() == 2) {
				for (int nx : pi.clauses[c].lit) {
					if (nx != x && pi.isLabel[litVariable(nx)] == VAR_UNDEFINED && BIGu[litNegation(nx)] != BIGIt) {
						rg[BIGid[x]].push_back(BIGid[litNegation(nx)]);
					}
				}
			}
		}
	}
	for (int x : cc) {
		BIGu[x] = BIGIt;
	}
	vector<int> ns;
	for (int i = 0; i < (int)cc.size(); i++) {
		BIGdfs2(g, i, ns);
	}
	int removed = 0;
	vector<vector<int> > sccs;
	vector<int> topo1;
	topo1.reserve(cc.size());
	for (int i = (int)ns.size() - 1; i >= 0; i--) {
		if (BIGu2[ns[i]] != 2) {
			vector<int> scc;
			BIGdfs3(rg, ns[i], scc);
			if (scc.size() > 1) {
				for (int& x : scc) {
					x = cc[x];
				}
				sccs.push_back(scc);
			}
			else {
				topo1.push_back(scc[0]);
			}
		}
	}
	for (auto& scc : sccs) {
		if (scc.size() > 1 && BIGu2[litVariable(scc[0])] != 3) {
			sort(scc.begin(), scc.end());
			bool unsat = false;
			for (unsigned j = 1; j < scc.size(); j++) {
				if (scc[j - 1] == litNegation(scc[j])) {
					unsat = true;
					break;
				}
			}
			if (unsat) {
				for (int l : scc) {
					vector<int> pCls = pi.litClauses[l];
					vector<int> nCls = pi.litClauses[litNegation(l)];
					for (int c : pCls) {
						pi.removeClause(c);
						rLog.removeClause(1);
					}
					for (int c : nCls) {
						pi.removeClause(c);
						rLog.removeClause(1);
					}
				}
				pi.addClause({scc[0]});
				pi.addClause({litNegation(scc[0])});
				rLog.removeVariable((int)scc.size() - 1);
				continue;
			}
			for (unsigned j = 0; j < scc.size(); j++) {
				assert(BIGu2[litVariable(scc[j])] != 3);
				BIGu2[litVariable(scc[j])] = 3;
			}
			for (unsigned j = 1; j < scc.size(); j++) {
				vector<int> pCls = pi.litClauses[scc[j]];
				vector<int> nCls = pi.litClauses[litNegation(scc[j])];
				for (int c : pCls) {
					pi.removeLiteralFromClause(scc[j], c);
				}
				for (int c : nCls) {
					pi.removeLiteralFromClause(litNegation(scc[j]), c);
				}
				for (int c : pCls) {
					bool f = false;
					bool fn = false;
					for (int l : pi.clauses[c].lit) {
						if (l == scc[0]) f = true;
						if (l == litNegation(scc[0])) fn = true;
					}
					if (fn) {
						pi.removeClause(c);
						rLog.removeClause(1);
					}
					else if (!f) {
						pi.addLiteralToClause(scc[0], c);
					}
				}
				for (int c : nCls) {
					bool f = false;
					bool fn = false;
					for (int l : pi.clauses[c].lit) {
						if (l == litNegation(scc[0])) f = true;
						if (l == scc[0]) fn = true;
					}
					if (fn) {
						pi.removeClause(c);
						rLog.removeClause(1);
					}
					else if (!f) {
						pi.addLiteralToClause(litNegation(scc[0]), c);
					}
				}
				trace.setEqual(scc[0], scc[j]);
			}
			rLog.removeVariable((int)scc.size() - 1);
			removed += (int)scc.size() - 1;
		}
	}
	if (sccs.size() > 0 || cc.size() == 1) {
		return removed;
	}
	if (!doTC) return 0;
	if (doneUnhiding) return 0;
	doneUnhiding = true;
	rLog.stopTechnique(Log::Technique::EE);
	rLog.startTechnique(Log::Technique::UH);
	int stamp = 1;
	vector<pair<int, int> > pcc(cc.size());
	for (int i = 0; i < (int)cc.size(); i++) {
		pcc[i] = {cc[i], i};
	}
	sort(pcc.begin(), pcc.end());
	
	vector<int> trueLits;
	vector<pair<int, int> > rmLit;
	vector<int> rmClause;
	vector<pair<int, int> > qrs;
	for (int i = 1; i < (int)pcc.size(); i++) {
		if (pcc[i].F == litNegation(pcc[i - 1].F)) {
			int a = pcc[i].S;
			int b = pcc[i - 1].S;
			qrs.push_back({a, b});
		}
	}
	vector<int> leIndex(cc.size());
	vector<int> riIndex(cc.size());
	vector<int> up1(cc.size());
	vector<int> up2(cc.size());
	vector<int> perm(cc.size());
	for (int i = 0; i < (int)cc.size(); i++) perm[i] = i;
	for (int trys = 0; trys < 9; trys++) { // magic constant
		random_shuffle(perm.begin(), perm.end());
		for (int i = 0; i < (int)cc.size(); i++) BIGu2[i] = 0;
		for (int i : perm) {
			genIndex(g, rg, i, i, i, stamp, leIndex, riIndex, up1, up2, trys%3);
		}
		// FAILED LITERALS
		for (auto qq : qrs) {
			if (BIGisPath(qq.F, qq.S, leIndex, riIndex, up1, up2)) trueLits.push_back(cc[qq.S]);
			else if (BIGisPath(qq.S, qq.F, leIndex, riIndex, up1, up2)) trueLits.push_back(cc[qq.F]);
		}
		for (int x : cc) {
			BIGu2[x] = 5;
		}
		//HLE
		for (int l1 : cc) {
			for (int c : pi.litClauses[l1]) {
// 				if (pi.clauses[c].lit.size() <= 2) continue;
				for (int l2 : pi.clauses[c].lit) {
					if (l2 != l1 && BIGu2[l2] == 5) {
						if (BIGisPath(BIGid[l1], BIGid[l2], leIndex, riIndex, up1, up2)) {
							rmLit.push_back({l1, c});
						}
					}
				}
			}
		}
		//HTE
		for (int l1 : cc) {
			for (int c : pi.litClauses[l1]) {
				if (pi.clauses[c].lit.size() <= 2) continue;
				for (int l2 : pi.clauses[c].lit) {
					if (l2 != l1 && BIGu2[litNegation(l2)] == 5) {
						if (BIGisPath(BIGid[litNegation(l2)], BIGid[l1], leIndex, riIndex, up1, up2)) {
							rmClause.push_back(c);
						}
					}
				}
			}
		}
		for (int x : cc) {
			BIGu2[x] = 0;
		}
	}
	sort(rmLit.begin(), rmLit.end());
	rmLit.erase(unique(rmLit.begin(), rmLit.end()), rmLit.end());
	for (auto rm : rmLit) {
		pi.removeLiteralFromClause(rm.F, rm.S);
		rLog.removeLiteral(1);
		removed++;
	}
	sort(rmClause.begin(), rmClause.end());
	rmClause.erase(unique(rmClause.begin(), rmClause.end()), rmClause.end());
	for (int c : rmClause) {
		pi.removeClause(c);
		rLog.removeClause(1);
	}
	sort(trueLits.begin(), trueLits.end());
	trueLits.erase(unique(trueLits.begin(), trueLits.end()), trueLits.end());
	for (int l : trueLits) {
		int rmC = 0;
		if (litValue(l) == true) {
			rmC = setVariable(litVariable(l), true);
		}
		else {
			rmC = setVariable(litVariable(l), false);
		}
		removed++;
		rLog.removeVariable(1);
		rLog.removeClause(rmC);
	}
	rLog.stopTechnique(Log::Technique::UH);
	rLog.startTechnique(Log::Technique::EE);
	return removed;
}

int Preprocessor::doBIG(bool doTC) {
	rLog.startTechnique(Log::Technique::EE);
	if (!rLog.requestTime(Log::Technique::EE)) {
		rLog.stopTechnique(Log::Technique::EE);
		return 0;
	}
	int removed = 0;
	BIGIt++;
	if ((int)BIGu.size() < 2*pi.vars) BIGu.resize(2*pi.vars);
	if ((int)BIGu2.size() < 2*pi.vars) BIGu2.resize(2*pi.vars);
	if ((int)BIGid.size() < 2*pi.vars) BIGid.resize(2*pi.vars);
	vector<int> checkLit = pi.tl.getBinaryLiterals("BIG");
	for (int lit : checkLit) {
		if (!rLog.requestTime(Log::Technique::EE)) break;
		if (BIGu[lit] == BIGIt || pi.isLabel[litVariable(lit)]) continue;
		removed += tryBIG(lit, doTC);
	}
	log(removed, " variables removed by EE");
	rLog.stopTechnique(Log::Technique::EE);
	return removed;
}

void Preprocessor::doBIG2(bool doTC) {
	BIGIt++;
	rLog.startTechnique(Log::Technique::EE);
	if ((int)BIGu.size() < 2*pi.vars) BIGu.resize(2*pi.vars);
	if ((int)BIGu2.size() < 2*pi.vars) BIGu2.resize(2*pi.vars);
	if ((int)BIGid.size() < 2*pi.vars) BIGid.resize(2*pi.vars);
	for (int lit = 0; lit < 2*pi.vars; lit++) {
		if (BIGu[lit] == BIGIt || pi.isLabel[litVariable(lit)]) continue;
		if (tryBIG(lit, doTC)) {
			print("fail BIG ", litToDimacs(lit));
			abort();
		}
	}
	rLog.stopTechnique(Log::Technique::EE);
}