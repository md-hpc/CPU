#include "avx.h"

fvec clipv(fvec v, float m) {
	__m256 min = _mm256_set1_ps(m);
	__m256 a = _mm256_load_ps(v);
	__m256 mask = _mm256_cmp_ps(a, min, _CMP_GT_OQ);

	a = _mm256_blendv_ps(a, min, mask);
	_mm256_store_ps(dst, a);
}

void sqrtv(fvec *dst, fvec *a) {
	__mm256 a;

	a = _mm256_load_ps(va);
	a = _mm256_sqrt_ps(a);
	_mm256_store_ps(dst,a);
}

fvec permute(fvec x) {
	pack p;
	for (int i = 0; i < VISZE; i++) {
		p.d[(i+1)%VSIZE] = x;
	}
	return p.v;
}
