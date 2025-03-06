#ifndef PARTICLE8_H
#define PARTICLE8_H

#include <vector>

#include "avx.h"
#include "particle.h"

using namespace std;

class vec8 {
public:
	vec8();
	vec8(f8 x, f8 y, f8 z);

	vec get(int i);
	void set(const vec &v, int i);

	vec8 operator+(const vec8 &other);
	
	vec8 operator*(const float c);
	vec8 operator*(const f8 c);
	
	vec8 &operator+=(const f8 c);
	vec8 &operator+=(const vec8 &other);

	vec8 &operator*=(const float c);
	vec8 &operator*=(const f8 c);
	
	void permutev();

	f8 x;
	f8 y;
	f8 z;
};

class particle8 {
public:
	particle8();
	particle8(const vec8 &r);
	particle8(const vec8 &r, const vec8 &v); 

	particle get(int i);
	void set(const particle &p, int i);

	vec8 r;
	vec8 v;
};

class p8buf {
public:
	p8buf();
	int append(const particle &sp);
	particle8 get();

	int i;
	particle8 p;
};

class particle8_vector {
public:
	void append(const particle &p);

	void set(const particle &p, int i); 
	particle get(int i); 

	void resize(int sz);
	int size1();
	int size8(); 
	
	
	particle8 &operator[](int i); 


private:
	vector<particle8> v;
	int n;
};



void printpv(vector<particle8_vector> &ps);
int pcount(vector<particle8_vector> &ps);
#endif
