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

	__m256 x = _mm256_set_ps(1,2,3,4,5,6,7,8);
	__m256 k = _mm256_set1_ps(2);
	x = _mm256_div_ps(x,k);

	__m256i n = _mm256_castps

	return 0;
}
