#include <iostream>
#include <sstream>
#include <vector>
#include <cstdint>

#include "outputreader.hpp"

using namespace std;
namespace maxPreprocessor {

int OutputReader::readSolution(istream& input) {
	status = 0;
	string line;
	while (getline(input, line)) {
		if (line[0] == 'o') {
			stringstream ss;
			ss<<line;
			string temp;
			ss>>temp;
			ss>>ansW;
		}
		if (line[0] == 's') {
			stringstream ss;
			ss<<line;
			string s;
			while (ss>>s) {
				if (s == "UNSATISFIABLE") {
					status = 2;
				}
				if (s == "OPTIMUM") {
					status = 1;
				}
			}
		}
		if (line[0] == 'v') {
			stringstream ss;
			ss<<line;
			string temp;
			ss>>temp;
			int lit;
			trueLits.clear();
			while (ss>>lit) {
				trueLits.push_back(lit);
			}
		}
	}
	if (status == 0) return 1;
	return 0;
}

}