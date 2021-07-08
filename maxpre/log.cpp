#include <vector>
#include <iostream>
#include <cassert>
#include <set>
#include <map>

#include "log.hpp"

using namespace std;
namespace maxPreprocessor {

LogT::LogT() {
	rClauses = 0;
	rVariables = 0;
	rLiterals = 0;
	rLabels = 0;
}

void LogT::print(ostream& out) {
	out<<"c "<<rClauses<<" "<<rVariables<<" "<<rLiterals<<" "<<rLabels<<'\n';
}

Log::Log() {
	activeTechnique = none;
	tLog.resize(30);
	tTimer.resize(30);
	timeLimit = 0;
	gatesExtracted = 0;
	labelsMatched = 0;
	binaryCoresFound = 0;
	toReallocate = 0;
	initialWeightRange = 0;
	weightRange = 0;
}

void Log::startTechnique(Technique t) {
	assert(activeTechnique == none);
	activeTechnique = t;
	tTimer[activeTechnique].start();
}

void Log::stopTechnique(Technique t) {
	assert(activeTechnique == t);
	tTimer[activeTechnique].stop();
	activeTechnique = none;
}

void Log::removeClause(int c) {
	if (c == 0) return;
	assert(activeTechnique != none);
	tLog[activeTechnique].rClauses+=c;
}

void Log::removeVariable(int c) {
	if (c == 0) return;
	assert(activeTechnique != none);
	tLog[activeTechnique].rVariables+=c;
}

void Log::removeLiteral(int c) {
	if (c == 0) return;
	assert(activeTechnique != none);
	tLog[activeTechnique].rLiterals+=c;
}

void Log::removeLabel(int c) {
	if (c == 0) return;
	assert(activeTechnique != none);
	tLog[activeTechnique].rLabels+=c;
}

void Log::print(ostream& out) {
	out<<"c removed clauses, variables, literals, labels\n";
	out<<"c BCE\n";
	tLog[BCE].print(out);
	out<<"c UP\n";
	tLog[UP].print(out);
	out<<"c BVE\n";
	tLog[BVE].print(out);
	out<<"c SE\n";
	tLog[SE].print(out);
	out<<"c SSR\n";
	tLog[SSR].print(out);
	out<<"c SLE\n";
	tLog[SLE].print(out);
	out<<"c BCR\n";
	tLog[BCR].print(out);
	out<<"c SIE\n";
	tLog[SIE].print(out);
	out<<"c EE\n";
	tLog[EE].print(out);
	out<<"c BVA\n";
	tLog[BVA].print(out);
	out<<"c GSLE\n";
	tLog[GSLE].print(out);
	out<<"c FLP\n";
	tLog[FLP].print(out);
	out<<"c UH\n";
	tLog[UH].print(out);
	out<<"c LS\n";
	tLog[LS].print(out);
	out.flush();
}

void Log::printTime(ostream& out) {
	out<<"c Preprocess time\n";
	out<<"c BCE "<<getTime(Technique::BCE)<<'\n';
	out<<"c UP "<<getTime(Technique::UP)<<'\n';
	out<<"c BVE "<<getTime(Technique::BVE)<<'\n';
	out<<"c SE "<<getTime(Technique::SE)<<'\n';
	out<<"c SSR "<<getTime(Technique::SSR)<<'\n';
	out<<"c SLE "<<getTime(Technique::SLE)<<'\n';
	out<<"c BCR "<<getTime(Technique::BCR)<<'\n';
	out<<"c SIE "<<getTime(Technique::SIE)<<'\n';
	out<<"c EE "<<getTime(Technique::EE)<<'\n';
	out<<"c BVA "<<getTime(Technique::BVA)<<'\n';
	out<<"c GSLE "<<getTime(Technique::GSLE)<<'\n';
	out<<"c FLP "<<getTime(Technique::FLP)<<'\n';
	out<<"c UH "<<getTime(Technique::UH)<<'\n';
	out<<"c LS "<<getTime(Technique::LS)<<'\n';
	out.flush();
}

void Log::printInfo(ostream& out) {
	out<<"c Labels matched "<<labelsMatched<<'\n';
	out<<"c Gates extracted "<<gatesExtracted<<'\n';
	out<<"c Binary cores found "<<binaryCoresFound<<'\n';
	out<<"c Original weight range "<<initialWeightRange<<'\n';
	out<<"c Preprocessed weight range "<<weightRange<<'\n';
}

double Log::getTime(Technique t) {
	return tTimer[t].getTime().count();
}

void Log::startTimer() {
	totalTimer.start();
}

void Log::stopTimer() {
	totalTimer.stop();
}

Log::Technique Log::charToTechnique(char c) {
	if (c == 'b') return Technique::BCE;
	else if (c == 'u') return Technique::UP;
	else if (c == 'v') return Technique::BVE;
	else if (c == 's') return Technique::SE;
	else if (c == 'r') return Technique::SSR;
	else if (c == 'l') return Technique::SLE;
	else if (c == 'c') return Technique::BCR;
	else if (c == 'i') return Technique::SIE;
	else if (c == 'e') return Technique::EE;
	else if (c == 'a') return Technique::BVA;
	else if (c == 'g') return Technique::GSLE;
	else if (c == 'p') return Technique::FLP;
	else if (c == 'h') return Technique::UH;
	else if (c == 't') return Technique::LS;
	else {
		assert(0);
	}
}

void Log::timePlan(double timeLimit_, string useTechniques) {
	eTime.resize(30);
	// magic constants
	eTime[Technique::BCE] = 0.05;
	eTime[Technique::UP] = 0.03;
	eTime[Technique::BVE] = 0.3;
	eTime[Technique::SE] = 0.2;
	eTime[Technique::SSR] = 0.15;
	eTime[Technique::SLE] = 0.03;
	eTime[Technique::BCR] = 0.03;
	eTime[Technique::SIE] = 0.05;
	eTime[Technique::EE] = 0.15;
	eTime[Technique::BVA] = 0.15;
	eTime[Technique::GSLE] = 0.08;
	eTime[Technique::FLP] = 0.15;
	eTime[Technique::UH] = 0.15;
	eTime[Technique::LS] = 0.15;
	set<Technique> useT;
	for (char c : useTechniques) {
		if (c == '[' || c == ']' || c == '#') continue;
		useT.insert(charToTechnique(c));
	}
	if (timeLimit_ < 0) timeLimit_ = 0;
	timeLimit += timeLimit_;
	double useTotal = 0;
	for (Technique t : useT) {
		assert(t < eTime.size());
		useTotal += eTime[t];
	}
	tTimeLimit.resize(eTime.size());
	for (Technique t : useT) {
		if (useTotal < 1e-9) tTimeLimit[t] = 0;
		else tTimeLimit[t] += timeLimit_*eTime[t]/useTotal;
	}
}

bool Log::requestTime(Technique t) {
	if (timeLimit > infTime/2) return true;
	assert(activeTechnique == t);
	if (askHistory.size() == 0 || askHistory.back() != t) newRequest(t);
	if (totalTimer.getTime().count() > timeLimit) return false;
	if (tTimer[t].getTime().count() < tTimeLimit[t]) return true;
	tTimeLimit[t] += toReallocate; // give reallocate time
	toReallocate = 0;
	return tTimer[t].getTime().count() < tTimeLimit[t];
}

void Log::newRequest(Technique t) {
	askHistory.push_back(t);
	if (totalTimer.getTime().count() > timeLimit) return; // timelimit anyway
	if (tTimer[t].getTime().count() < tTimeLimit[t]) return; // has time left
	// try to allocate more time
	map<Technique, double> excessTime;
	for (int i = (int)askHistory.size() - 2; i >= 0; i--) {
		if (askHistory[i] == t) break;
		if (excessTime.count(askHistory[i]) == 0) {
			Technique tt = askHistory[i];
			if (tTimeLimit[tt] - tTimer[tt].getTime().count() > 0) {
				excessTime[tt] = tTimeLimit[tt] - tTimer[tt].getTime().count();
			}
		}
	}
	for (auto et : excessTime) {
		// magic constant
		tTimeLimit[t] += et.second/3; // reallocate third of the excess
		tTimeLimit[et.first] -= et.second/3;
	}
}

void Log::neverAgain(char c) {
	Technique t = charToTechnique(c);
	toReallocate += max((double)0, tTimeLimit[t] - tTimer[t].getTime().count());
}

bool Log::isTimeLimit() {
	return timeLimit < infTime/2;
}

}