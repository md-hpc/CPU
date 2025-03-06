#include <stdio.h>

typedef void (*fptr)(int);

void p(int a) {
	printf("%d\n",a);
}

int main() {
	fptr ptr = p;

	p(5);
}
