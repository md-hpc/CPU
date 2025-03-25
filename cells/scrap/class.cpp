#include <stdio.h>

class A { 
public:
	A(int c) {
		if (c == 1) {
			a = 4;
		} else {
			a = 3;
		}
	}

	const int a;
};


int main(int argc, char **argv) {
	A obj = A(atoi(argv[1]));

	printf("%d\n",obj.a);
}
