#ifndef MAXPP_PREPROCESSOR_HPP
#define MAXPP_PREPROCESSOR_HPP

#include <vector>
#include <queue>
#include <set>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "global.hpp"
#include "preprocessedinstance.hpp"
#include "probleminstance.hpp"
#include "trace.hpp"
#include "timer.hpp"
#include "clause.hpp"
#include "log.hpp"
#include "AMSLEX.hpp"

namespace maxPreprocessor {
// Variables are indexed 0..n-1
// Positive literal is v*2, negative is v*2+1
class Preprocessor {
public:
	ProblemInstance pi;
	int originalClauses;
	int originalVars;
	
	Log rLog;
	
	int logLevel = 1;
	bool printComments = true;
	
	Trace trace;
	
	Preprocessor(const std::vector<std::vector<int> >& clauses_, const std::vector<uint64_t>& weights_, uint64_t topWeight_);
	
	bool isTautology(const Clause& clause) const;
	
	// Returns number of clauses removed
	int setVariable(int var, bool value);
	
	// This is called only in the beginning since no tautologies are added
	void removeTautologies();
	
	int eliminateReduntantLabels();
	// This is called only in the beginning
	void identifyLabels();
	
	// This is called only in the beginning
	void createLabels();
	
	const uint64_t polyHashMul = 1000000007; // magic constant
	
	int tryUP(int lit);
	int tryUPAll();
	
	int doUP();
	void doUP2();
	
	int removeDuplicateClauses();
	
	// Is clause a subsumed by clause b?
	// Supposes that a and b are sorted
	bool isSubsumed(const std::vector<int>& a, const std::vector<int>& b) const;
	
	AMSLEX amsLex;
	void trySEHash(std::vector<int>& clauses, int tLit, std::vector<int>& toRemove);
	void trySEAmsLex(std::vector<int>& clauses, std::vector<int>& toRemove);
	int trySESlow(int lit);
	void trySE(std::vector<int>& clauses, std::vector<int>& toRemove);
	void trySEgen(int lit, std::vector<int>& toRemove);
	
	int doSE();
	void doSE2();
	
	bool BVEgate;
	void printC(int c) const;
	std::pair<std::vector<int>, int> searchXor(int var) const;
	std::pair<std::vector<int>, int> searchITE(int var) const;
	std::pair<std::vector<int>, int> searchAndOr(int dLit, int defVars, const std::set<int>& binaryClauseLits) const;
	// Dont give labels to this
	int tryBVE(int var);
	int tryBVE2(int var);
	std::vector<int> tryBVEGE(int var);
	
	std::vector<uint64_t> getBVEHash(const std::vector<int>& cs, int var, int sw) const;
	
	int doBVE();
	void doBVE2();
	
	
	// Supposes that the clauses are in sorted order
	// Can we do SSR such that var is removed from c2?
	bool canSSR(int var, const Clause& c1, const Clause& c2);
	
	
	bool SSRC(int c1, int c2, int var);
	int trySSRAmsLex(int var);
	int trySSRHash(int var);
	int trySSR2(int var);
	int trySSR(int var);
	int trySSRgen(int var);
	
	int doSSR();
	void doSSR2();
	
	int tryBCE(int lit);
	
	int doBCE();
	void doBCE2();
	
	bool vSubsumed(std::vector<int>& v1, std::vector<int>& v2);
	int trySLESlow(int lb1, int lb2);
	int doSLE();
	void doSLE2();
	
	bool tryBCR(int c, int l11);
	int doBCR();
	void doBCR2();
	
	bool SIErndCheck(int litX, int litY);
	int try2SIE(int litX, int litY);
	int trySIE(int lit);
	int doSIE();
	void doSIE2();
	
	void genIndex(std::vector<std::vector<int> >& g, std::vector<std::vector<int> >& rg, int x, int u1, int u2, int& stamp, std::vector<int>& le, std::vector<int>& ri, std::vector<int>& up1, std::vector<int>& up2, int order);
	int BIGIt;
	std::vector<int> BIGu, BIGu2, BIGid;
	void BIGdfs1(int x, std::vector<int>& ns);
	void BIGdfs2(std::vector<std::vector<int> >& g, int x, std::vector<int>& ns);
	void BIGdfs3(std::vector<std::vector<int> >& rg, int x, std::vector<int>& scc);
	bool BIGisPath(int x, int to, std::vector<int>& leIndex, std::vector<int>& riIndex, std::vector<int>& up1, std::vector<int>& up2);
	int tryBIG(int lit, bool doTC);
	int doBIG(bool doTC);
	void doBIG2(bool doTC);
	bool doneUnhiding;
	
	std::vector<uint64_t> sfH;
	std::vector<uint64_t> tMul;
	std::unordered_map<uint64_t, int> BVAHashTable;
	void addBVAHash(std::vector<int>& lits, std::unordered_map<uint64_t, int>& hashes);
	int canBVA(int c, int d, int lit);
	int tryBVA(int lit, std::unordered_map<uint64_t, int>& hashes);
	int doBVA();
	void doBVA2();
	
	void GSLEBT(int i, uint64_t w, std::vector<int>& sel, std::vector<uint64_t>& weights, std::vector<std::vector<int> >& hs, bool& found, uint64_t& itLim);
	bool GSLEtryBackTrack(std::vector<std::vector<int> >& hs, std::vector<uint64_t>& weights, uint64_t w, uint64_t itLim);
	int tryGSLE(int lb);
	int doGSLE();
	void doGSLE2();
	
	int tryFLP(std::vector<int> fLit, int clause);
	int doFLP();
	
	void tryLFF(int lb);
	void findLabeledFormula();
	
	int tryLS(int lbl);
	void tryLSBCE(int lit, std::unordered_set<int>& deletedClauses, std::unordered_set<int>& touchedList, std::vector<std::pair<int, int> >& blockedClauses);
	int doLS();
	
	void CBIGdfs1(int x, std::vector<int>& ns, std::vector<std::pair<int, std::pair<int, int> > >& condEdges);
	int findConditionalGraph(int lit, std::vector<std::pair<int, std::pair<int, int> > >& condEdges);
	std::vector<std::pair<int, std::pair<int, int> > > findConditionalComponents();
	
	int removeEmptyClauses();
	
	bool validTechniques(std::string techniques) const;
	bool validPreTechniques(std::string techniques) const;
	
	PreprocessedInstance getPreprocessedInstance();
	
	int doPreprocess(const std::string& techniques, int l, int r, bool debug, bool topLevel);
	void preprocess(std::string techniques, double timeLimit, bool debug, bool BVEgate, bool initialCall, bool matchLabels);
	
	int skipTechnique;
	
	std::mt19937 randGen;
	
	template<typename T>
	void log(T t) {
		if (logLevel == 0) return;
		std::cerr << t << std::endl;
	}

	template<typename T, typename... Args>
	void log(T t, Args... args) {
		if (logLevel == 0) return;
		std::cerr << t;
		log(args...);
	}
	
	template<typename T>
	void print(T t) {
		if (!printComments) return;
		std::cout << t << std::endl;
	}

	template<typename T, typename... Args>
	void print(T t, Args... args) {
		if (!printComments) return;
		std::cout << t;
		print(args...);
	}
	
	template<typename T>
	T getRand(T lo, T hi) {
		return std::uniform_int_distribution<T>(lo, hi)(randGen);
	}
};
}
#endif