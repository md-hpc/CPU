class vec {
	vec(float x, float y, float z) x(x), y(y), z(z) {}
	vec() : x(NAN), y(NAN), z(NAN) {}

	vec operator+(const vec &other) {
		return vec(x + other.x, y + other.y, z + other.z);
	}

	vec operator*(float k) {
		return vec(k * x, k * y, k * z);
	}

	float x;
	float y;
	float z;
}
