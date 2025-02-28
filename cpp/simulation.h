#ifndef I_SIMULATION
#define I_SIMULATION

#include <vector>
#include "avx.h"
#include "fvv.h"
#include "common.h"

using namespace std;

class voxel {
public:
	voxel(int i, int j, int k);

	int i;
	int j;
	int k;
};

typedef struct {
	void (*worker)(int);
	int core;
	int start;
	int stop;
} worker_spec_t;

class simulation {
public:
	simulation(int argc, char **argv);

	void simulate();

	fvec ljv(fvec x, fvec y, fvec z);
	fvec submv(fvec va, fvec vb);
	ivec cellv(fvec x, fvec y, fvec z);
	fvec apbcfv(fvec x);
	
	void velocity_update_worker(int hci);
	void position_update_worker(int hci);
	void cell_update_worker(int hci);
	
	voxel voxelof(int i);
	int cell(int i, int j, int k);
	
	void *run_worker(void *arg);
	void thread(void (*worker)(void*), int threads);

	vector<vector<particle8>> cells;
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
	int RESOLUTION;
	int THREADS;

	int LOGFD;
	int FD;

	int t;

	char default_log[32] = "validate/cells-vec";
	char default_path[32] = "viz/particles";

	// "constants" for LJ computation
	float EP4;
	float SPSS;
	float TPST;

};


#endif
