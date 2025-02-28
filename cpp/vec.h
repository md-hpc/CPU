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
