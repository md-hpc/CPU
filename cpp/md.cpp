#include <stdio.h>
#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

using namespace std;

ALGO = ALGO_CELLS_VEC;

typedef int ivec __attribute__ ((vector_size(32)));
typedef float fvec __attribute__((vector_size(32)));
#define VBYTES sizeof(fvec)
#define VSIZE (sizeof(fvec)/sizeof(float))

typedef union {
	fvec v;
	float d[VSIZE];
} pack;

typedef union {
	ivec v;
	int d[VSIZE];
} ipack;


#define VAI(v,i) (((pack*)v)->d[i])
#define VI(v,i) (((pack*)&(v))->d[i])
#define VK(x) {k,k,k,k,k,k,k,k}


// performs no bounds checking.
class fvec_vector_vec {
public:
	fvec_vector_vec() : vi(0) {}

	vec get(int i) {
		return vec(
			VAI(x,i),
			VAI(y,i),
			VAI(z,i)
		);
	}

	void set(vec v, int i) {
		VAI(x,i) = v.x;
		VAI(y,i) = v.y;
		VAI(z,i) = v.z;
	}

	void append(vec v) {
		const float zero = 0;	
		if (n % 8 == 0) {
			x.push_back(zero);
			y.push_back(zero);
			z.push_back(zero);
		}
		set(v,n++);
	}

	void resize(int sz) {
		if (sz > n) {
			prinf("Tried to increase size of fvec_vector_vec")
			throw 1;
		}

		n = sz;
		r = VSIZE - n % VSIZE;
		for (int i = n; i < n + r; i++) {
			set(vec(0,0,0),i);
		}
	}
	
	int size() {
		return n;
	}

	vector<fvec> x;
	vector<fvec> y;
	vector<fvec> z;

private:
	int n;
}



simulation s;

int main(int argc, char **argv) {
	simulation s;
	int np;
	particle *p;

	parse_cli(argc, argv);

	init_particles(particles);

	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				if (di==0 && dj==0 && dk==0)
					continue;
				offsets.push_back((voxel_index(di,dj,dk));
			}
		}
	}

	np = particles.size();
	for (int pi = 0; pi < np; pi++) {
		p = &particles[pi];
		s.r[p->cell].append(p->r);
		s.v[p->cell].append(p->v);
	}

	for (t = 0; t < TIMESTEPS) {
		thread(velocity_update, N_CELL);
		thread(motion_update, N_CELL);
		thread(cell_update, N_CELL);
	}
}

void permute(fvec *x, int n) {
	pack tmp, *dst;
	tmp.v = *x;
	dst = (pack*) x;
	for (int i = 0; i < VISZE; i++) {
		dst->d[(i+n)%VSIZE] = tmp.d[i];
	}
}



