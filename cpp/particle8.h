class vec8 {
	vec8();
	vec8(fvec x, fvec y, fvec z);

	vec get(int i);
	void set(vec v, int i);

	fvec x;
	fvec y;
	fvec z;
}

class particle8 {
	particle8();
	particle8(vec8 r);
	particle8(vec8 r, vec8 v); 

	particle get(int i):
	void set(particle p, int i);

	vec8 r;
	vec8 p;
}

class particle8_vector {
public:
	void append(const &particle p);

	void set(const &particle p, int i); 
	particle get(int i); 

	particle8 &operator[](int i); 

private:
	vector<particle8> v;
	int n;
}
