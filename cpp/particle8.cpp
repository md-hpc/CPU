#include <math.h>

vec8::vec8() {
	const float nan = NAN;
	x = nan;
	y = nan;
	z = nan;
}

vec8::vec8(fvec x, fvec y, fvec z) : x(x), y(y), z(z) {}

vec vec8::get(int i) {
	return vec(
		((pack)x).d[i],
		((pack)y).d[i],
		((pack)z).d[i]
	)
}

vec vec8::set(const vec &v, int i) {
	((pack)x).d[i] = v.x;
	((pack)y).d[i] = v.y;
	((pack)z).d[i] = v.z;
}

particle8::particle8() {
	r = vec8();
	v = vec8();
}

particle8::particle8(const vec8 &r) : r(r) {
	v = vec8();
}

particle8::particle8(const vec8 &r, const vec8 &v) : r(r), v(v) {}

particle particle8::get(int i) {
	return particle(
		r.get(i),
		v.get(i)
	);
}

particle particle8::set(const particle &p, int i) {
	r.set(p.r,i);
	v.set(p.v,i);
}

void particle8_vector::append(const &particle p) {
	if (n % VSIZE == 0) {
		v.push_back(particle8());
	}
	set(p,n);
	n++;
}

void particle8_vector::set(const &particle p, int i) {
	v[n/VSIZE].set(p, n%VSIZE);
}

particle particle8_vector::get(int i) {
	return v[n/VSIZE].get(i%VSIZE);
}

particle8 &particle8_vector::operator[](int i) {
	return v[i];
}
