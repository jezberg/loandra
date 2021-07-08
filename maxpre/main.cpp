#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>
#include <sstream>
#include <map>

#include "preprocessorinterface.hpp"
#include "inputreader.hpp"
#include "outputreader.hpp"
#include "timer.hpp"
#include "utility.hpp"
using namespace std;

map<string, string> getFlags(int argc, char* argv[]) {
	map<string, string> ret;
	for (int i = 3; i < argc; i++) {
		string s(argv[i]);
		if (s.size() == 0 || s[0] != '-') {
			cout<<"c Invalid arg "<<s<<endl;
			cerr<<"Invalid arg "<<s<<endl;
		}
		else {
			int p = -1;
			int prd = 0;
			for (int j = 0; j < (int)s.size(); j++) {
				if (s[j] == '=') {
					p = j;
					break;
				}
				if (j == prd && s[j] == '-') prd++;
			}
			if (p == -1) {
				ret[s.substr(prd)] = "";
			}
			else {
				ret[s.substr(prd, p-prd)] = s.substr(p+1);
			}
		}
	}
	return ret;
}

string parseTechniques(map<string, string>& flags) {
	string techniques = "[bu]#[buvsrgc]";
	if (flags.count("techniques")) {
		techniques = flags["techniques"];
		cout<<"c Techniques "<<techniques<<endl;
		cerr<<"Techniques "<<techniques<<endl;
	}
	else {
		cout<<"c No -techniques= given, defaulting to "<<techniques<<endl;
		cerr<<"No -techniques= given, defaulting to "<<techniques<<endl;
	}
	return techniques;
}

double parseTimeLimit(map<string, string>& flags) {
	double timeLimit = 1e9;
	if (flags.count("timelimit")) {
		stringstream ss;
		ss<<flags["timelimit"];
		ss>>timeLimit;
		cout<<"c Preprocess time limit "<<timeLimit<<endl;
		cerr<<"Preprocess time limit "<<timeLimit<<endl;
	}
	else {
		cout<<"c No -timelimit= given, defaulting to inf"<<endl;
		cerr<<"No -timelimit= given, defaulting to inf"<<endl;
	}
	return timeLimit;
}

bool parseBVEgate(map<string, string>& flags) {
	bool BVEgate = false;
	if (flags.count("bvegate")) {
		if (flags["bvegate"] == "1") {
			BVEgate = true;
			cout<<"c BVE gate extraction enabled"<<endl;
			cerr<<"BVE gate extraction enabled"<<endl;
		}
		else if (flags["bvegate"] == "0") {
			BVEgate = false;
			cout<<"c BVE gate extraction disabled"<<endl;
			cerr<<"BVE gate extraction disabled"<<endl;
		}
		else {
			cout<<"Invalid bvegate flag"<<endl;
			cerr<<"Invalid bvegate flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -bvegate= given, defaulting to disabled"<<endl;
		cerr<<"No -bvegate= given, defaulting to disabled"<<endl;
	}
	return BVEgate;
}

bool parseLabelMatching(map<string, string>& flags) {
	bool labelMatching = false;
	if (flags.count("matchlabels")) {
		if (flags["matchlabels"] == "1") {
			labelMatching = true;
			cout<<"c Label matching enabled"<<endl;
			cerr<<"Label matching enabled"<<endl;
		}
		else if (flags["matchlabels"] == "0") {
			labelMatching = false;
			cout<<"c Label matching disabled"<<endl;
			cerr<<"Label matching disabled"<<endl;
		}
		else {
			cout<<"Invalid matchlabels flag"<<endl;
			cerr<<"Invalid matchlabels flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -matchlabels given, defaulting to disabled"<<endl;
		cerr<<"No -matchlabels given, defaulting to disabled"<<endl;
	}
	return labelMatching;
}

bool parseProblemType(map<string, string>& flags) {
	bool maxSat = true;
	if (flags.count("problemtype")) {
		for (char& c : flags["problemtype"]) {
			c = tolower(c);
		}
		if (flags["problemtype"] == "sat") {
			maxSat = false;
			cout<<"c Problem type is SAT"<<endl;
			cerr<<"Problem type is SAT"<<endl;
		}
		else if (flags["problemtype"] == "maxsat" || flags["problemtype"] == "max-sat") {
			maxSat = true;
			cout<<"c Problem type is Max-SAT"<<endl;
			cerr<<"Problem type is Max-SAT"<<endl;
		}
		else {
			cout<<"Invalid problemtype flag"<<endl;
			cerr<<"Invalid problemtype flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -problemtype given, defaulting to Max-SAT"<<endl;
		cerr<<"No -problemtype given, defaulting to Max-SAT"<<endl;
	}
	return maxSat;
}

int parseSkipTechnique(map<string, string>& flags) {
	int skipTechnique = 0;
	if (flags.count("skiptechnique")) {
		stringstream ss;
		ss<<flags["skiptechnique"];
		ss>>skipTechnique;
		cout<<"c Skiptechnique "<<skipTechnique<<endl;
		cerr<<"Skiptechnique "<<skipTechnique<<endl;
		
		if (skipTechnique <= 0 || skipTechnique > 1000000000) {
			cout<<"Invalid skiptechnique flag"<<endl;
			cerr<<"Invalid skiptechinque flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -skiptechnique given, defaulting to disabled"<<endl;
		cerr<<"No -skiptechnique given, defaulting to disabled"<<endl;
	}
	return skipTechnique;
}

int parseVerb(map<string, string>& flags) {
	int verb = 1;
	if (flags.count("verb")) {
		if (flags["verb"] == "0") {
			verb = 0;
		}
		else if (flags["verb"] == "1") {
			verb = 1;
		}
		else {
			cout<<"Invalid verb flag"<<endl;
			cerr<<"Invalid verb flag"<<endl;
			exit(0);
		}
		cout<<"c Verb "<<verb<<endl;
		cerr<<"Verb "<<verb<<endl;
	}
	else {
		cout<<"c No -verb given, defaulting to 1"<<endl;
		cerr<<"No -verb given, defaulting to 1"<<endl;
	}
	return verb;
}

void printHelp(ostream& out) {
	out<<endl;
	out<<"The first argument is the instance file, the second is preprocess, reconstruct or solve."<<endl;
	out<<endl;
	
	out<<"An example of using the preprocessor:"<<endl;
	out<<"\t./preprocessor test.wcnf preprocess -techniques=[bu]#[buvsrg] -mapfile=test.map > preprocessed.wcnf"<<endl;
	out<<"\t./solver < preprocessed.wcnf > sol0.sol"<<endl;
	out<<"\t./preprocessor sol0.sol reconstruct -mapfile=test.map > solution.sol"<<endl;
	out<<endl;
	out<<"Another way to do the same thing:"<<endl;
	out<<"\t./preprocessor test.wcnf solve -solver=./solver -techniques=[bu]#[buvsrg] > solution.sol"<<endl;
	out<<endl;
	
	out<<"-techniques (default: [bu]#[buvsrgc])"<<endl;
	out<<"\tstring:"<<endl;
	out<<"\tThis string defines the preprocessing techniques to use and the order of them."<<endl;
	out<<"\tEach letter corresponds to a preprocessing technique. Each preprocessing technique is applied until its fixpoint."<<endl;
	out<<"\tTechniques inside brackets are applied until all of them are in fixpoint. The brackets work recursively. "<<endl;
	out<<"\tIf # character is given, all techniques before it are applied before group detection and adding labels (techniques available before labeling are BCE, UP and SE)."<<endl;
	out<<"\tTechniques:"<<endl;
	out<<"\tb = blocked clause elimintation"<<endl;
	out<<"\tu = unit propagation"<<endl;
	out<<"\tv = bounded variable elimination"<<endl;
	out<<"\ts = subsumption elimination"<<endl;
	out<<"\tr = self subsuming resolution"<<endl;
	out<<"\tl = subsumed label elimintion"<<endl;
	out<<"\tc = binary core removal"<<endl;
	out<<"\ta = bounded variable addition"<<endl;
	out<<"\tg = generalized subsumed label elimination"<<endl;
	out<<"\te = equivalence elimination"<<endl;
	out<<"\th = unhiding techniques (failed literals, hidden tautology elimination, hidden literal elimination)"<<endl;
	out<<"\tt = structure labeling"<<endl;
	out<<"\tp = failed label probing"<<endl;
	out<<endl;
	
	out<<"-solver (default: disabled)"<<endl;
	out<<"\tstring:"<<endl;
	out<<"\tThe solver to use to solve the preprocessed instance"<<endl;
	out<<endl;
	
	out<<"-solverflags (default: disabled)"<<endl;
	out<<"\tstring:"<<endl;
	out<<"\tThe flags to use with the solver"<<endl;
	out<<"\tFor example -solver=./LMHS -solverflags=\"--infile-assumps --no-preprocess\" results in using the command ./LMHS preprocessed.wcnf --infile-assumps --no-preprocess > sol0.sol"<<endl;
	out<<endl;
	
	out<<"-mapfile (default: disabled)"<<endl;
	out<<"\tstring:"<<endl;
	out<<"\tThe file to write the solution reconstruction map"<<endl;
	out<<endl;
	
	out<<"-problemtype (default: maxsat)"<<endl;
	out<<"\tstring: {maxsat, sat}"<<endl;
	out<<"\tShould the problem be preprocessed as a MaxSAT or SAT instance"<<endl;
	out<<endl;
	
	out<<"-outputformat (default: wpms)"<<endl;
	out<<"\tstring: {original}"<<endl;
	out<<"\tBy default the preprocessor always gives the output in weighted partial MaxSAT format"<<endl;
	out<<"\tOutput in SAT format by setting this to original when preprocessing SAT instances"<<endl;
	out<<endl;
	
	out<<"-timelimit (default: inf)"<<endl;
	out<<"\tdouble: [0, 500000000]"<<endl;
	out<<"\tLimit for preprocessing time in seconds"<<endl;
	out<<endl;
	
	out<<"-skiptechnique (default: disabled)"<<endl;
	out<<"\tint: [1, 1000000000]"<<endl;
	out<<"\tSkip a preprocessing technique if it seems to be not effective in x tries (x is given in this flag)"<<endl;
	out<<"\tRecommended values for this could be something between 10 and 1000"<<endl;
	out<<endl;
	
	out<<"-matchlabels (default: 0)"<<endl;
	out<<"\tbool: {0, 1}"<<endl;
	out<<"\tUse label matching technique to reduce the number of labels"<<endl;
	out<<endl;
	
	out<<"-bvegate (default: 0)"<<endl;
	out<<"\tbool: {0, 1}"<<endl;
	out<<"\tUse BVE gate extraction to extend BVE"<<endl;
	out<<"\tNote: applying BCE will destroy all recognizable gate structures"<<endl;
	out<<endl;
	
	out<<"-verb (default: 1)"<<endl;
	out<<"\tint: [0, 1]"<<endl;
	out<<"\tIf verb is 0 the preprocessor will output less stuff to the standard error"<<endl;
	out<<endl;
}

int isHelp(char* arg) {
	return string(arg) == "-help" || string(arg) == "--help" || string(arg) == "-h" || string(arg) == "--h";
}

int main(int argc, char* argv[]){
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	if ((argc > 1 && isHelp(argv[1])) || (argc > 2 && isHelp(argv[2])) || (argc > 3 && isHelp(argv[3]))) {
		printHelp(cout);
		return 0;
	}
	if (argc < 3) {
		cout<<"Use -help"<<endl;
		cout<<"Use -help"<<endl;
		return 0;
	}
	auto flags = getFlags(argc, argv);
	if (flags.count("h") || flags.count("help")) {
		printHelp(cout);
	}
	string type(argv[2]);
	assert(type == "solve" || type == "preprocess" || type == "reconstruct");
	string file(argv[1]);
	if (type == "reconstruct") {
		if (!flags.count("mapfile")) {
			cout<<"Give mapfile with -mapfile= flag"<<endl;
			cerr<<"Give mapfile with -mapfile= flag"<<endl;
			return 0;
		}
		string mapFile = flags["mapfile"];
		maxPreprocessor::OutputReader opr;
		ifstream in(file);
		int readStatus = opr.readSolution(in);
		in.close();
		if (readStatus > 0) {
			cout<<"Failed to parse solution"<<endl;
			cerr<<"Failed to parse solution"<<endl;
			return 0;
		}
		if (opr.status == 2) {
			cout<<"s UNSATISFIABLE"<<endl;
		}
		else {
			ifstream mapF(mapFile);
			int vars, ppVars, originalVars;
			mapF>>vars>>ppVars>>originalVars;
			vector<int> varMap(vars);
			for (int i = 0; i < vars; i++) mapF>>varMap[i];
			vector<int> trueLits;
			for (int lit : opr.trueLits) {
				if (abs(lit) > vars) continue;
				if (lit > 0) lit = varMap[abs(lit)-1];
				else lit = -varMap[abs(lit)-1];
				trueLits.push_back(lit);
			}
			maxPreprocessor::Trace trace;
			int traceLines;
			mapF>>traceLines;
			trace.operations.resize(traceLines);
			trace.data.resize(traceLines);
			for (int i = 0; i < traceLines; i++) {
				mapF>>trace.operations[i];
				int sz;
				mapF>>sz;
				trace.data[i].resize(sz);
				for (int j = 0; j < sz; j++) {
					mapF>>trace.data[i][j];
				}
			}
			trace.printSolution(cout, trueLits, opr.ansW, ppVars, originalVars);
		}
		return 0;
	}
	
	if (type == "solve") {
		// just check that -solver flag is used
		if (!(flags.count("solver") && flags["solver"].size() > 0)) {
			cout<<"Please specify the solver"<<endl;
			cerr<<"Please specify the solver"<<endl;
			return 0;
		}
	}
	
	string techniques = parseTechniques(flags);
	double timeLimit = parseTimeLimit(flags);
	bool BVEgate = parseBVEgate(flags);
	bool labelMatching = parseLabelMatching(flags);
	bool maxSat = parseProblemType(flags);
	int skipTechnique = parseSkipTechnique(flags);
	
	ifstream instanceFile(file);
	if (instanceFile.fail()) {
		cout<<"Failed to read the input file"<<endl;
		cerr<<"Failed to read the input file"<<endl;
		return 0;
	}
	maxPreprocessor::InputReader inputReader;
	int readStatus = inputReader.readClauses(instanceFile, maxSat);
	instanceFile.close();
	
	if (readStatus > 0) {
		cout<<"Failed to parse input instance: "<<inputReader.readError<<endl;
		cerr<<"Failed to parse input instance: "<<inputReader.readError<<endl;
		return 0;
	}
	
	int outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS;
	if (flags["outputformat"] == "original") {
		outputFormat = inputReader.inputFormat;
		if (outputFormat == maxPreprocessor::INPUT_FORMAT_MS) {
			// preprocessor works in labeled cnf so it cannot output pure maxsat
			outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS;
		}
		string outf;
		if (outputFormat == maxPreprocessor::INPUT_FORMAT_WPMS) {
			outf = "weighted partial Max-SAT";
		}
		else if (outputFormat == maxPreprocessor::INPUT_FORMAT_SAT) {
			outf = "SAT";
		}
		else {
			return 0;
		}
		cout<<"c Outputformat "<<outf<<endl;
		cerr<<"Outputformat "<<outf<<endl;
	}
	
	int verb = parseVerb(flags);
	
	maxPreprocessor::PreprocessorInterface pif(inputReader.clauses, inputReader.weights, inputReader.top);
	pif.setBVEGateExtraction(BVEgate);
	pif.setLabelMatching(labelMatching);
	pif.setSkipTechnique(skipTechnique);
	
	maxPreprocessor::Timer preprocessTimer;
	preprocessTimer.start();
	
	pif.preprocess(techniques, verb, timeLimit);
	
	preprocessTimer.stop();
	
	cerr<<"Preprocess time: "<<preprocessTimer.getTime().count()<<endl;
	if (verb > 0) pif.printTimeLog(cerr);
	if (type == "preprocess") {
		if (verb > 0) pif.printTechniqueLog(cerr);
		pif.printTechniqueLog(cout);
		pif.printTimeLog(cout);
		pif.printInfoLog(cout);
		pif.printInstance(cout, outputFormat);
		
		if (flags.count("mapfile")) {
			string mapFile = flags["mapfile"];
			cout<<"c Outputting reconstruction map to "<<mapFile<<endl;
			cerr<<"Outputting reconstruction map to "<<mapFile<<endl;
			ofstream out(mapFile);
			pif.printMap(out);
			out.close();
		}
		else {
			cout<<"c No -mapfile= given, will not ouput reconstruction map"<<endl;
			cerr<<"No -mapfile= given, will not ouput reconstruction map"<<endl;
		}
	}
	if (type == "solve") {
		string solver;
		if (flags.count("solver") && flags["solver"].size() > 0) {
			solver = flags["solver"];
			cout<<"c Using solver "<<solver<<endl;
			cerr<<"Using solver "<<solver<<endl;
			cout<<"c Solver flags "<<flags["solverflags"]<<endl;
			cerr<<"Solver flags "<<flags["solverflags"]<<endl;
		}
		else {
			cout<<"Please specify the solver"<<endl;
			cerr<<"Please specify the solver"<<endl;
			return 0;
		}
		
		ofstream out("preprocessed.wcnf");
		pif.printInstance(out, outputFormat);
		out.close();
		
		maxPreprocessor::Timer solveTimer;
		solveTimer.start();
		if (system((solver + " preprocessed.wcnf " + flags["solverflags"] + " > sol0.sol").c_str())) {
			cout<<"Solver error"<<endl;
			cerr<<"Solver error"<<endl;
			return 0;
		}
		solveTimer.stop();
		
		maxPreprocessor::OutputReader opr;
		ifstream in("sol0.sol");
		readStatus = opr.readSolution(in);
		in.close();
		if (readStatus > 0) {
			cout<<"Failed to parse solution"<<endl;
			cerr<<"Failed to parse solution"<<endl;
			return 0;
		}
		
		if (opr.status == 2) {
			cout<<"s UNSATISFIABLE"<<endl;
		}
		else {
			pif.printSolution(opr.trueLits, cout, opr.ansW);
		}
		cerr<<"Preprocess time: "<<preprocessTimer.getTime().count()<<", Solve time: "<<solveTimer.getTime().count()<<endl;
		if (verb > 0) pif.printTimeLog(cerr);
	}
}