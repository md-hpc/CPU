#ifndef I_FFV
#define I_FFV 1 

#include <vector>
#include "avx.h"
#include "common.h"

using namespace std;

// performs no bounds checking.
class fvec_vector_vec {
public:
	fvec_vector_vec() : n(0) {}

	vec get(int i);
	void set(const vec &v, int i); 

	vec8 get8(int i);
	void set(const vec8 &v, int i);

	void append(vec v);
	void resize(int sz);	
	int size(); 

	vector<fvec> x;
	vector<fvec> y;
	vector<fvec> z;

private:
	int n;
};

#endif
