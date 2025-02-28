#include <vector>

// performs no bounds checking.
class fvec_vector_vec {
public:
	fvec_vector_vec() : n(0) {}

	vec get(int i);
	void set(vec v, int i); 

	void append(vec v);
	void resize(int sz);	
	int size(); 

	vector<fvec> x;
	vector<fvec> y;
	vector<fvec> z;

private:
	int n;
}


