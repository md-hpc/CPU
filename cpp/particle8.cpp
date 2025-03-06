#include <math.h>
#include <stdio.h>

#include "particle8.h"

vec8::vec8() {
	f8 nan = VK(NAN);
	x = nan;
	y = nan;
	z = nan;
}

vec8::vec8(f8 x, f8 y, f8 z) : x(x), y(y), z(z) {}

vec vec8::get(int i) {
	return vec(
		VI(x,i),
		VI(y,i),
		VI(z,i)
	);
}

void vec8::set(const vec &v, int i) {
	VI(x,i) = v.x;
	VI(y,i) = v.y;
	VI(z,i) = v.z;
}

vec8 vec8::operator+(const vec8 &other) {
    return vec8(x + other.x, y + other.y, z + other.z);
}

vec8 vec8::operator*(const f8 c) {
    return vec8(c*x, c*y, c*z);
}

vec8 vec8::operator*(const float c) {
	return vec8(c*x, c*y, c*z);
}

vec8 &vec8::operator*=(const f8 c) {
    x *= c;
    y *= c;
    z *= c;
   
    return *this;
}

vec8 &vec8::operator*=(const float c) {
	x *= c;
	y *= c;
	z *= c;

	return *this;
}

vec8 &vec8::operator+=(const f8 c) {
    x += c;
    y += c;
    z += c;

    return *this;
}

vec8 &vec8::operator+=(const vec8 &other) {
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;
}

void vec8::permutev() {
	x = permute(x);
	y = permute(y);
	z = permute(z);
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
		v.get(i),
		-1
	);
}

void particle8::set(const particle &p, int i) {
	r.set(p.r,i);
	v.set(p.v,i);
}

p8buf::p8buf() : i(0) {}

int p8buf::append(const particle &sp) {
	if (i == 8) {
		printf("p8buf overflow");
		throw 1;
	}
	p.set(sp,i++);
	return i == 8;
}

particle8 p8buf::get() {
	particle8 ret = p;
	i = 0;
	p = particle8();
	return ret;
}

void particle8_vector::append(const particle &p) {
	if (n % VSIZE == 0) {
		v.push_back(particle8());
	}
	set(p,n);
	n++;
}

void particle8_vector::set(const particle &p, int i) {
	v[n/VSIZE].set(p, n%VSIZE);
}

particle particle8_vector::get(int i) {
	return v[n/VSIZE].get(i%VSIZE);
}

particle8 &particle8_vector::operator[](int i) {
	return v[i];
}

void particle8_vector::resize(int sz) {
	particle p; // NAN initalized
	n = sz;
	if (n % 8) {
		for (int i = n; i < n + (8 - n % 8); i++) {
			set(p,i);
		}
	}
	v.resize(sz / 8 + (sz % 8 != 0));
}

int particle8_vector::size1() {
	return n;
}

int particle8_vector::size8() {
	return v.size();
}

void printpv(vector<particle8_vector> &ps) {
	int n = ps.size();
	for (int i = 0; i < n; i++) {
		int np = ps[i].size1();
		if (np == 0)
			continue;
		printf("\t%d {", i);
		for (int j = 0; j < np; j++) {
			particle p = ps[i].get(j);
			printf("(%.1f %.1f %.1f, %.1e %.1e %.1e) ", p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
		}
		printf("}\n");
	}
}
