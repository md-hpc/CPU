class simulation {
public:
	int parse_cli(int argc, char **argv);

	fvec ljv(fvec r);
	fvec submv(fvec va, fvec vb);
	fvec cellv(fvec x, fvec y, fvec z);
	fvec apbcfv(fvec x);
	
	void velocity_update_worker(int hci, int tid);
	void position_update_worker(int hci, int tid);
	void cell_update_worker(int hci, int tid);

	int voxelof(int i, int j, int k);

	void *run_worker(void *arg);
	void thread(void (*worker)(void*), int threads);

	vector<fvec_vector_vec> r;
	vector<fvec_vector_vec> v;
	vector<vector<vec>> vos;
	vector<vector<vec>> ros;
	vector<vector<int>> ocis;

	float SIGMA;
	float EPSILON;
	float DT;
	int UNIVERSE_SIZE;
	int PARTICLES;
	float CUTOFF;

	float LJ_MIN;
	float L;
	int N_CELL;

	int TIMESTEPS;
	int SEED;
	int RESOLUTION;
	int THREADS;
};

class voxel {
	voxel(int i, int j, int k);

	int i;
	int j;
	int k;
}

typedef struct {
	void (*worker)(int);
	int core;
	int start;
	int stop;
} worker_spec_t;
