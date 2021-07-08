#include <cstdint>
#include <vector>
#include <algorithm>
#include <cassert>

#include "preprocessor.hpp"
#include "preprocessedinstance.hpp"
#include "preprocessorinterface.hpp"
#include "global.hpp"
#include "utility.hpp"

#define F first
#define S second

using namespace std;
namespace maxPreprocessor {
	PreprocessorInterface::PreprocessorInterface(const vector<vector<int> >& clauses, const vector<uint64_t>& weights, uint64_t topWeight_, bool inProcessMode_) 
	: preprocessor(clauses, weights, topWeight_), topWeight(topWeight_), inProcessMode(inProcessMode_) {
		variables = 0;
		originalVariables = 0;
		for (auto& clause : clauses) {
			for (int lit : clause) {
				assert(lit != 0);
				variables = max(variables, abs(lit));
			}
		}
		originalVariables = variables;
		preprocessed = false;
		useBVEGateExtraction = false;
		useLabelMatching = false;
		skipTechnique = 0;
		if (inProcessMode) {
			PPVarToSolverVar.resize(variables);
			solverVarToPPVar.resize(variables);
			for (auto& clause : clauses) {
				for (int lit : clause) {
					PPVarToSolverVar[abs(lit)-1] = abs(lit);
					solverVarToPPVar[abs(lit)-1] = abs(lit);
				}
			}
		}
	}
	
	uint64_t PreprocessorInterface::getRemovedWeight() {
		return preprocessor.trace.removedWeight;
	}
	
	void PreprocessorInterface::preprocess(string techniques, int logLevel, double timeLimit) {
		preprocessor.logLevel = logLevel;
		preprocessor.printComments = false;
		preprocessor.skipTechnique = skipTechnique;
		
		preprocessor.preprocess(techniques, timeLimit, false, useBVEGateExtraction, !preprocessed, useLabelMatching);
		preprocessed = true;
		variables = max(variables, preprocessor.pi.vars);
	}
	
	int PreprocessorInterface::getNumClauses() {
		return preprocessor.getPreprocessedInstance().clauses.size();
	}
	
	void PreprocessorInterface::getInstance(std::vector<std::vector<int> >& retClauses, std::vector<uint64_t>& retWeights, std::vector<int>& retLabels) {
		preprocessedInstance = preprocessor.getPreprocessedInstance();
		
		retClauses = preprocessedInstance.clauses;
		retWeights = preprocessedInstance.weights;
		
		for (unsigned i = 0; i < preprocessedInstance.labels.size(); i++) {
			retLabels.push_back(litToSolver(litToDimacs(preprocessedInstance.labels[i].F)));
		}
		for (auto& clause : retClauses) {
			for (int& lit : clause) {
				lit = litToDimacs(lit);
				variables = max(variables, abs(lit));
				lit = litToSolver(lit);
			}
		}
		
		for (uint64_t& w : retWeights) {
			if (w == HARDWEIGHT) {
				w = topWeight;
			}
		}
	}
	
	vector<int> PreprocessorInterface::reconstruct(const vector<int>& trueLiterals) {
		vector<int> ppTrueLiterals;
		for (int lit : trueLiterals) {
			lit = litToPP(lit);
			if (lit == 0) continue;
			ppTrueLiterals.push_back(lit);
		}
		return preprocessor.trace.getSolution(ppTrueLiterals, 0, variables, originalVariables).F;
	}
	
	void PreprocessorInterface::printSolution(const vector<int>& trueLiterals, ostream& output, uint64_t ansWeight) {
		vector<int> ppTrueLiterals;
		for (int lit : trueLiterals) {
			lit = litToPP(lit);
			if (lit == 0) continue;
			ppTrueLiterals.push_back(lit);
		}
		preprocessor.trace.printSolution(output, ppTrueLiterals, ansWeight, variables, originalVariables);
	}
	
	uint64_t PreprocessorInterface::getTopWeight() {
		return topWeight;
	}
	
	
	
	int PreprocessorInterface::addVar(int var) {
		if (!inProcessMode) return 0;
		assert(var >= 0);
		if (var == 0) {
			var = solverVarToPPVar.size() + 1;
		}
		if (var > (int)solverVarToPPVar.size()) {
			solverVarToPPVar.resize(var);
		}
		if (solverVarToPPVar[var-1] != 0) return 0;
		int nv = preprocessor.pi.addVar()+1;
		solverVarToPPVar[var-1] = nv;
		assert((int)PPVarToSolverVar.size() < nv);
		PPVarToSolverVar.resize(nv);
		PPVarToSolverVar[nv-1] = var;
		return var;
	}
	
	int PreprocessorInterface::addLabel(int lbl, uint64_t weight) {
		if (!inProcessMode) return 0;
		lbl = addVar(lbl);
		if (lbl == 0) return 0;
		if (weight >= topWeight) weight = HARDWEIGHT;
		int iVar = solverVarToPPVar[lbl-1]-1;
		preprocessor.pi.addClause({negLit(iVar)}, weight);
		preprocessor.pi.isLabel[iVar] = VAR_FALSE;
		return lbl;
	}
	
	bool PreprocessorInterface::alterWeight(int lbl, uint64_t weight) {
		if (!inProcessMode) return false;
		if (lbl < 1) return false;
		int iVar = solverVarToPPVar[lbl-1]-1;
		assert(iVar >= 0);
		int softClause = -1;
		if (preprocessor.pi.isLabel[iVar] == VAR_TRUE) {
			assert(preprocessor.pi.litClauses[posLit(iVar)].size() == 1);
			softClause = preprocessor.pi.litClauses[posLit(iVar)][0];
		}
		else if (preprocessor.pi.isLabel[iVar] == VAR_FALSE) {
			assert(preprocessor.pi.litClauses[negLit(iVar)].size() == 1);
			softClause = preprocessor.pi.litClauses[negLit(iVar)][0];
		}
		else {
			return false;
		}
		if (weight >= topWeight) {
			preprocessor.pi.isLabel[iVar] = 0;
			weight = HARDWEIGHT;
		}
		assert(!preprocessor.pi.clauses[softClause].isHard());
		preprocessor.pi.clauses[softClause].weight = weight;
		if (weight == HARDWEIGHT) {
			assert(preprocessor.pi.clauses[softClause].isHard());
		}
		else {
			assert(preprocessor.pi.labelWeight(iVar) == weight);
		}
		return true;
	}
	
	bool PreprocessorInterface::addClause(std::vector<int> clause) {
		if (!inProcessMode) return false;
		for (int& lit : clause) {
			assert(lit != 0);
			int tlit = litToPP(lit);
			if (tlit == 0) {
				int nvar = addVar(abs(lit));
				assert(nvar == abs(lit));
				tlit = litToPP(lit);
				assert(lit == litToSolver(tlit));
			}
			tlit = litFromDimacs(tlit);
			lit = tlit;
		}
		sort(clause.begin(), clause.end());
		clause.erase(unique(clause.begin(), clause.end()), clause.end());
		for (int i = 1; i < (int)clause.size(); i++) {
			if (clause[i] == litNegation(clause[i-1])) return true;
		}
		preprocessor.pi.addClause(clause);
		return true;
	}
	
	bool PreprocessorInterface::labelToVar(int lbl) {
		if (!inProcessMode) return false;
		if (lbl < 1) return false;
		int iVar = solverVarToPPVar[lbl-1]-1;
		assert(iVar >= 0);
		int softClause = -1;
		if (preprocessor.pi.isLabel[iVar] == VAR_TRUE) {
			assert(preprocessor.pi.litClauses[posLit(iVar)].size() >= 1);
			softClause = preprocessor.pi.litClauses[posLit(iVar)][0];
		}
		else if (preprocessor.pi.isLabel[iVar] == VAR_FALSE) {
			assert(preprocessor.pi.litClauses[negLit(iVar)].size() >= 1);
			softClause = preprocessor.pi.litClauses[negLit(iVar)][0];
		}
		else {
			return true;
		}
		preprocessor.pi.isLabel[iVar] = 0;
		preprocessor.pi.removeClause(softClause);
		return true;
	}
	
	bool PreprocessorInterface::resetRemovedWeight() {
		if (!inProcessMode) return false;
		preprocessor.trace.removedWeight = 0;
		return true;
	}
	
	void PreprocessorInterface::setBVEGateExtraction(bool use) {
		useBVEGateExtraction = use;
	}
	
	void PreprocessorInterface::setLabelMatching(bool use) {
		useLabelMatching = use;
	}
	
	void PreprocessorInterface::setSkipTechnique(int value) {
		skipTechnique = value;
	}
	
	int PreprocessorInterface::litToSolver(int lit) {
		if (PPVarToSolverVar.size() < abs(lit)) PPVarToSolverVar.resize(abs(lit));
		if (PPVarToSolverVar[abs(lit)-1] == 0) {
			PPVarToSolverVar[abs(lit)-1] = solverVarToPPVar.size() + 1;
			solverVarToPPVar.push_back(abs(lit));
		}
		if (lit > 0) return PPVarToSolverVar[abs(lit)-1];
		else return -PPVarToSolverVar[abs(lit)-1];
	}
	
	int PreprocessorInterface::litToPP(int lit) {
		if ((int)solverVarToPPVar.size() <= abs(lit)-1) return 0;
		if (lit > 0) return solverVarToPPVar[abs(lit)-1];
		else return -solverVarToPPVar[abs(lit)-1];
	}
	
	void PreprocessorInterface::printInstance(ostream& output, int outputFormat) {
		std::vector<std::vector<int> > clauses;
		std::vector<uint64_t> weights;
		std::vector<int> labels;
		getInstance(clauses, weights, labels);
		
		assert(outputFormat == INPUT_FORMAT_WPMS || outputFormat == INPUT_FORMAT_SAT);
		
		if (clauses.size() == 0) {
			clauses.push_back({-1, 1});
			weights.push_back(topWeight);
		}
		
		if (outputFormat == INPUT_FORMAT_WPMS) {
			output<<"c assumptions ";
			for (unsigned i = 0; i < labels.size(); i++) {
				output<<labels[i];
				if (i + 1 < labels.size()) {
					output<<" ";
				}
			}
			output<<"\n";
		}
		
		if (outputFormat == INPUT_FORMAT_WPMS) {
			output<<"p wcnf "<<max((int)solverVarToPPVar.size(), 1)<<" "<<clauses.size()<<" "<<topWeight<< '\n';
			for (unsigned i = 0; i < clauses.size(); i++) {
				output<<weights[i]<<" ";
				for (int lit : clauses[i]) {
					output<<lit<<" ";
				}
				output<<"0\n";
			}
		}
		else if (outputFormat == INPUT_FORMAT_SAT) {
			output<<"p cnf "<<max((int)solverVarToPPVar.size(), 1)<<" "<<clauses.size()<<'\n';
			for (unsigned i = 0; i < clauses.size(); i++) {
				for (int lit : clauses[i]) {
					output<<lit<<" ";
				}
				output<<"0\n";
			}
		}
		else {
			abort();
		}
		output.flush();
	}
	void PreprocessorInterface::printTechniqueLog(ostream& output) {
		preprocessor.rLog.print(output);
	}
	void PreprocessorInterface::printTimeLog(ostream& output) {
		preprocessor.rLog.printTime(output);
	}
	void PreprocessorInterface::printInfoLog(ostream& output) {
		preprocessor.rLog.printInfo(output);
	}
	void PreprocessorInterface::printMap(ostream& output) {
		output<<solverVarToPPVar.size()<<" "<<variables<<" "<<originalVariables<<'\n';
		for (int t : solverVarToPPVar) {
			output<<t<<" ";
		}
		output<<'\n';
		output<<preprocessor.trace.operations.size()<<'\n';
		for (unsigned i = 0; i < preprocessor.trace.operations.size(); i++) {
			output<<preprocessor.trace.operations[i]<<" ";
			output<<preprocessor.trace.data[i].size()<<" ";
			for (int a : preprocessor.trace.data[i]) {
				output<<a<<" ";
			}
			output<<'\n';
		}
		output.flush();
	}
}


