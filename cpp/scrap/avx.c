#include <immintrin.h>
#include <stdio.h>
#include <math.h>

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
		printf("%f ", d[i]);
	}
	printf("\n");

	__m256 x = _mm256_set_ps(1,NAN,1,NAN,1,1,1,1);
	__m256 z = _mm256_set1_ps(0);
	__m256 nanm = _mm256_cmp_ps(x,x,_CMP_UNORD_Q);
	x = _mm256_blendv_ps(x,z,nanm);
	
	d = (float*) &x;
	for (int i = 0; i < 8; i++) {
		printf("%f ", d[i]);
	}
	printf("\n");

	return 0;
}
