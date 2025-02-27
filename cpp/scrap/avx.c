#include <immintrin.h>
#include <stdio.h>


typedef float fvec __attribute__ ((vector_size(32)));

int main(int argc, char **argv) {
	float k = atof(argv[1]);
	
	fvec a = {1,2,3,4,5,6,7,8};

	float *d = (float*) &a;

	for (int i = 0; i < 8; i++) {
		printf("%f ",d[i]);
	}
	printf("\n");
	return 0;
}
