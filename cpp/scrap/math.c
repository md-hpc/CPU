#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline float apbcf(float x, float m) {
	x = fmodf(x,m);
	return x < 0 ? x + m : x;
}

int main() {
	printf("%f\n",apbcf(-3,2.5));
}
