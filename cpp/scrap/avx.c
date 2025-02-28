#include <immintrin.h>
#include <stdio.h>


typedef float fvec __attribute__ ((vector_size(32)));

#define VK(x) {x,x,x,x,x,x,x,x}


fvec add(fvec a, fvec b) {
	return a + b;
}

int main(int argc, char **argv) {
	float k;

	if (argc != 2) {
		k = 1.;
	} else {
		k = atof(argv[1]);
	}
	const float K = k; 

	fvec a = {1,2,3,4,5,6,7,8};
	a += K;

	float *d = (float*) &a;

	for (int i = 0; i < 8; i++) {
		printf("%f ",d[i]);
	}
	printf("\n");
	return 0;
}
