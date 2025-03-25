#include <vector>
#include <unordered_map>
#include <cstdio>

using namespace std;

typedef unordered_map<int, int> dict;

int main() {
	vector<dict> d;

	d.resize(5);

	d[0][1] = 5;
	d[1][3] = 5;

	return 0;
}
