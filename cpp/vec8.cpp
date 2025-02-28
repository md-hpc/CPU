#include "vec8.h"

vec8::vec8(fvec x, fvec y, fvec z) : x(x), y(y), z(z) {}

vec8::vec8() {
	const float zero = 0;
	x = 0;
	y = 0;
	z = 0;
}

vec8 vec8::operator+(const vec8 &other) {
    return vec8(x + other.x, y + other.y, z + other.z);
}

vec8 vec8::operator*(const fvec c) {
    return vec8(c*x,c*y,c*z);
}

vec8 &vec8::operator*=(const fvec c) {
    x *= c;
    y *= c;
    z *= c;
   
    return *this;
}

vec8 &vec8::operator+=(const fvec c) {
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

vec8 &vec8::operator=(const vec8& other) {
    x = other.x;
    y = other.y;
    z = other.z;

    return *this;
}


