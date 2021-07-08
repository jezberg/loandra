#include <vector>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "clause.hpp"
#include "global.hpp"

using namespace std;
namespace maxPreprocessor {

Clause::Clause (vector<int> literals_, uint64_t weight_) : lit(literals_), weight(weight_) {
	updateHash();
}

bool Clause::isHard() const {
	return weight == HARDWEIGHT;
}

void Clause::updateHash() {
	hash = 0;
	for (uint64_t l : lit) {
		hash |= ((uint64_t)1 << (l&(uint64_t)63));
	}
}

void Clause::addLiteral(int l) {
	lit.push_back(0);
	bool f = false;
	for (int i = 0; i < (int)lit.size() - 1; i++) {
		if (lit[i] > l) {
			for (int j = (int)lit.size() - 1; j > i; j--) {
				lit[j] = lit[j - 1];
			}
			lit[i] = l;
			f = true;
			break;
		}
	}
	if (!f) lit.back() = l;
	
	updateHash();
}

void Clause::removeLiteral(int l) {
	bool f = false;
	
	for (int i = 0; i < (int)lit.size(); i++) {
		if (lit[i] == l) {
			f = true;
			for (int j = i + 1; j < (int)lit.size(); j++) {
				lit[j - 1] = lit[j];
			}
			lit.pop_back();
			break;
		}
	}
	
	assert(f == true);
	updateHash();
}

}