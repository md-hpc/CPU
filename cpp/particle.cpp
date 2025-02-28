#include <math.h>
#include "particle.h"

vec::vec() {
	x = nan;
	y = nan;
	z = nan;
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
