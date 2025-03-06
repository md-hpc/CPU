#include <math.h>
#include "particle.h"

vec::vec() {
	x = NAN;
	y = NAN;
	z = NAN;
}

vec::vec(float x, float y, float z) : x(x), y(y), z(z) {}

particle::particle() {
	r = vec();
	v = vec();
	cell = -1;
}

particle::particle(vec r) : r(r) {
	v = vec();
	cell = -1;
}

particle::particle(vec r, vec v, int cell) : r(r), v(v), cell(cell) {}
