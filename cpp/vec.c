vec fvec_vector_vec::get(int i) {
	return vec(
		VAI(x,i),
		VAI(y,i),
		VAI(z,i)
	);
}

void fvec_vector_vec::set(vec v, int i) {
	VAI(x,i) = v.x;
	VAI(y,i) = v.y;
	VAI(z,i) = v.z;
}

void fvec_vector_vec::append(vec v) {
	const float zero = 0;	
	if (n % 8 == 0) {
		x.push_back(zero);
		y.push_back(zero);
		z.push_back(zero);
	}
	set(v,n++);
}

void fvec_vector_vec::resize(int sz) {
	if (sz > n) {
		printf("Tried to increase size of fvec_vector_vec")
		throw 1;
	}

	n = sz;
	r = VSIZE - n % VSIZE;
	for (int i = n; i < n + r; i++) {
		set(vec(0,0,0),i);
	}
}

int fvec_vector_vec::size() {
	return n;
}

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
