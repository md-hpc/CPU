#ifdef PARTICLE_H
#define PARTICLE_H

class vec {
	vec();
	vec(float x, float y, float z);

	float x;
	float y;
	float z;
};

class particle {
	particle();
	partcile(vec r);
	particle(vec r, vec v, int cell);

	vec r;
	vec v;
	int cell;
};

#endif
