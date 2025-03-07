#ifndef I_SIMULATION
#define I_SIMULATION

#include <vector>
#include <pthread.h>

#include "avx.h"
#include "particle8.h"
#include "common.h"

using namespace std;


class simulation;

typedef struct {
	int core;
	int start;
	int stop;
	simulation *s;
} worker_spec_t;

class voxel {
public:
	voxel(int i, int j, int k);

	int i;
	int j;
	int k;
};

class simulation {
public:
	simulation(int argc, char **argv);

	void simulate();
	void save();

	vec8 lj(const vec8 &r, const vec8 &n);
	vec8 submv(const vec8 &va, const vec8 &vb);
	f8 submv(f8 va, f8 vb);
	i8 cellv(const vec8 &r);
	void apbcfv(vec8 &r);
	f8 apbcfv(f8 x);
	int apbci(int i);

	void velocity_update_worker(int hci, worker_spec_t* spec);
	void position_update_worker(int hci, worker_spec_t* spec);
	void cell_update_worker(int hci, worker_spec_t* spec);
	
	voxel voxelof(int i);
	int cell(int i, int j, int k);
	int cell(const vec &v);


	static void *run_worker(void *arg);
	void do_work(worker_spec_t *spec);
	void create_workers();
	void join_workers();

	vector<particle8_vector> particles;
	vector<vector<particle>> outbounds;
	vector<vector<int>> cios;

	float SIGMA;
	float EPSILON;
	float DT;
	int UNIVERSE_SIZE;
	int PARTICLES;
	float CUTOFF;

	float LJ_MIN;
	float L;
	int CELLS;

	int TIMESTEPS;
	int SEED;
	int THREADS;

	int LOGFD;
	int FD;
	int RESOLUTION;
	int SAVE;

	int t;

	char default_log[32] = "validate/cells-vec";
	char default_path[32] = "viz/particles";

	// "constants" for LJ computation
	float EP4;
	float SPSS;
	float TPST;

	pthread_t *tids;
	pthread_barrier_t barrier;
	pthread_barrier_t parent_barrier;
	worker_spec_t *specs;
};

#endif
