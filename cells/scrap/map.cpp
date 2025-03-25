#include <unordered_map>
#include <cstdio>

using namespace std;

typedef unordered_map<int, int> idict;

class A {
public:
	A() : a(1), b(2) {}

	int a;
	int b;
};

int main() {

	idict dict;

	dict[1] = 11;
	dict[3] = 13;
	dict[5] = 15;

	for (idict::iterator it = dict.begin(); it != dict.end(); it++) {
		printf("%d\n",*it);
	}
	return 0;

}
