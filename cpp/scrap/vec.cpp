#include <vector>
#include <stdio.h>

typedef float fvec __attribute__((vector_size(32)));

using namespace std;

int main() {
	fvec v = {1,2,3,4,5,6,7,8};

	vector<fvec> a;

	a.push_back(v);

	printf("%p\n",&a[0]);
	printf("%d\n",((long)&a[0])%32);
	printf("%d\n",sizeof(fvec));
}
