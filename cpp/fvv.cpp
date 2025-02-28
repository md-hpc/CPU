vec fvec_vector_vec::get(int i) {
	return vec(
		VAI(x,i),
		VAI(y,i),
		VAI(z,i)
	);
}

void fvec_vector_vec::set(const vec &v, int i) {
	VAI(x,i) = v.x;
	VAI(y,i) = v.y;
	VAI(z,i) = v.z;
}

vec8 fvec_vector_vec::get8(int i) {
	return vec8(
		x[i],
		y[i],
		z[i]
	);
}

void fvec_vector_vec::set8(const vec8 &v, int i) {
	x[i] = v.x;
	y[i] = v.y;
	z[i] = v.z;
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


