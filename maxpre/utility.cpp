#include "utility.hpp"
#include "global.hpp"

using namespace std;
namespace maxPreprocessor {

int litFromDimacs(int lit) {
	if (lit < 0) {
		return (-lit)*2-1;
	}
	else if(lit > 0) {
		return lit*2-2;
	}
	else {
		return -1;
	}
}

int litToDimacs(int lit) {
	if (lit&1) {
		return -(lit/2+1);
	}
	else {
		return lit/2+1;
	}
}

}