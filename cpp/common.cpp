#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "common.h"
#include <unistd.h>

using namespace std;

// Constants for liquid argon
/*
float SIGMA = 0.34;
float EPSILON = 120;
float CUTOFF = .85;
float M = 39.984 * 1e-3 / (6.022e23);
float DT = 1e-15;
*/
// Uhh
float SIGMA = 1;
float EPSILON = 1;
float CUTOFF = 2.5;
float M = 1;
float DT = 1e-4;

// length of a cell. Is CUTOFF if ALGO_CELLS and 1.2*CUTOFF if ALGO_LISTS
float R;


int UNIVERSE_SIZE = 5;
int N_PARTICLE = -1;
int N_TIMESTEP = 10;
int SEED = 0;
int RESOLUTION = 10;
int NEIGHBOR_REFRESH_RATE = 18;
int BR = -1;
int BN = -1;
int THREADS = 128;
int SAVE = 0;
mdalgo_t ALGO = ALGO_NONE;

float LJ_MIN; // minimum lj. Arbitrarily 4x the magnitude of the potential well
float CSQ; // cutoff squared

vec::vec(float x, float y, float z) : x(x), y(y), z(z) {};
vec::vec() : x(0), y(0), z(0) {};

vec vec::operator+(const vec &other) {
    return vec(x + other.x, y + other.y, z + other.z);
}

vec vec::operator*(const float c) {
    return vec(c*x,c*y,c*z);
}

vec &vec::operator*=(const float c) {
    x *= c;
    y *= c;
    z *= c;
    
    return *this;
}

vec &vec::operator+=(const float c) {
    x += c;
    y += c;
    z += c;

    return *this;
}

vec &vec::operator+=(const vec &other) {
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;


vec &vec::operator=(const vec& other) {
    x = other.x;
    y = other.y;
    z = other.z;

    return *this;
}

vec vec::operator%(const vec &other) {
    return vec(
        subm(other.x, x),
        subm(other.y, y),
        subm(other.z, z)
    );
}

static inline float apbcf(float x) {
    x = fmodf(x,L);
    return x < 0 ? x + L : x;
}

void vec::apbc() {
	x = apbcf(x);
	y = apbcf(y);
	z = apbcf(z);
}

int vec::cell() {
    return linear_idx(
        (int)x/R,
        (int)y/R,
        (int)z/R
    );
}

float vec::normsq() {
    return x*x+y*y+z*z;
}

float vec::norm() {
    return sqrt(x*x+y*y+z*z);
}

void vec::read(float *buf) {
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
}

#ifdef DEBUG
void vec::print() {
    printf("%.3f %.3f %.3f\n",x,y,z);
}

char *vec::str() {
    sprintf(strbuf,"(%.3f %.3f %.3f)",x,y,z);
    return strbuf;
}
#endif

int particle::counter = 0;

particle::particle() : r(vec(0,0,0)), v(vec(0,0,0)) {}
particle::particle(vec r) : r(r), v(vec(0,0,0)) {
    cell = r.cell();
	id = counter++;
}

int particle::interact(particle *pn) {
	float f, r;
	vec dv;

	if (pn == this) {
		return 0;
	}

	dv = pn->r % this->r;

	#ifdef DEBUG
	// printf("%s - %s = %s\n",pn->r.str(), this->r.str(), dv.str());
	#endif

	
	if (ALGO != ALGO_LISTS && dv.x < 0)
		return 0;

	r = dv.normsq();

	if (r > CSQ) 
		return 0;
	
	r = sqrt(r);
	
	f = LJ(r);
	if (f < -4*LJ(R_MAX)) {
		f = -4*LJ(R_MAX);
	}

	dv *= f/r * DT;
	
	v += dv;
	dv *= -1;
	pn->v += dv;
	
	return 1;
}


void particle::update_position() {
	vec dr = v * DT;

	#ifdef DEBUG
	if (dr.x > CUTOFF || dr.y > CUTOFF || dr.z > CUTOFF) {
		printf("%d has large velocity of %e", id, dr.norm());
	}
	#endif

	r += dr;
	cell = r.cell();
}

#ifdef DEBUG
char *particle::str() {
    sprintf(dbstr,"%d %d %d", id, cell, old_cell);
    return dbstr;
}
#endif



timer::timer() : time(0), running(false) {}

void timer::start() {
    last = rdtsc();
}

void timer::stop() {
    time += rdtsc() - last;
}

unsigned long timer::get() {
    return time;
}

static inline int apbci(int i) {
    i = i < 0 ? i + UNIVERSE_SIZE : i;
    i = i >= UNIVERSE_SIZE ? i - UNIVERSE_SIZE : i;
    return i;
}

int linear_idx(int i, int j, int k) {
    i = apbci(i);
    j = apbci(j);
    k = apbci(k);
    return i + j*UNIVERSE_SIZE + k*UNIVERSE_SIZE*UNIVERSE_SIZE;
}

void cubic_idx(int *res, int idx) {
    res[2] = idx / (UNIVERSE_SIZE * UNIVERSE_SIZE);
    idx -= res[2] * UNIVERSE_SIZE * UNIVERSE_SIZE;
    res[1] = idx / UNIVERSE_SIZE;
    idx -= res[1] * UNIVERSE_SIZE;
    res[0] = idx;
}

float subm(float a, float b) {
    // distance from a to b under periodic boundary condition of length L
    float opts[3] = {
        b - (a - L),
        b - a,
        b - (a + L)
    };
    float ds[3] = {
        abs(opts[0]),
        abs(opts[1]),
        abs(opts[2])
    };
    int mi = 0;
    mi = ds[mi] > ds[1] ? 1 : mi;
    mi = ds[mi] > ds[2] ? 2 : mi;
    return opts[mi];
}

float lj(float r) {
   	float f = LJ(r);
	if (f < -4*LJ(R_MAX)) {
		f = -4*LJ(R_MAX);
	}

	return f;
}

float frand() {
    return ((float)rand()/(float)RAND_MAX);
}


typedef struct {
    int core;
    int start;
    int stop;
    void (*kernel)(int, int);
} job_t;

void *run_kernel(void *spec) {
    job_t *job = (job_t*) spec;

    cpu_set_t cpuset;

    CPU_ZERO(&cpuset);
    CPU_SET(job->core, &cpuset);
    if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
        perror("could not set affinity");
        return (void*) 1;
    }
    
    for (int i = job->start; i < job->stop; i++) {
        job->kernel(i, job->core);
    }
    
    return (void*) 0;
}

void thread(void (*kernel)(int, int), int n) {
    pthread_t tids[THREADS];
    job_t jobs[THREADS];
    void *ret;
    int t;
    int bsize = (n + THREADS - 1) / THREADS;
    
    for (t = 0; t < THREADS; t++) {
        jobs[t].core = t;
        jobs[t].start = t * bsize;
        jobs[t].stop = (t + 1) * bsize < n ? (t + 1) * bsize : n;
        jobs[t].kernel = kernel;
        
        if (pthread_create(&tids[t], NULL, run_kernel, &jobs[t])) {
			perror("Could not start thread");
			exit(1);
		}
    }

    for (t = 0; t < THREADS; t++) {
        pthread_join(tids[t], &ret);
        if (ret) {
            exit(1);
        }
    }
}

int parse_cli(int argc, char **argv) {    
    if (ALGO == ALGO_NONE) {
		printf("Must set ALGO before calling parse_cli");
		return 1;
	}

	for (char **arg = &argv[1]; arg < &argv[argc]; arg+=2) {
        if (!strcmp(arg[0],"--sigma")) {
            SIGMA = atof(arg[1]);
        } else if (!strcmp(arg[0],"--epsilon")) {
            EPSILON = atof(arg[1]);
        } else if (!strcmp(arg[0],"--cutoff")) {
            CUTOFF = atof(arg[1]);
        } else if (!strcmp(arg[0],"--universe-size")) {
            UNIVERSE_SIZE = atol(arg[1]);
        } else if (!strcmp(arg[0],"--particles")) {
            N_PARTICLE = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--timesteps")) {
            N_TIMESTEP = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--dt")) {
            DT = atof(arg[1]);
        } else if (!strcmp(arg[0],"--seed")) {
            SEED = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--resolution")) {
            RESOLUTION = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--br")) {
            BR = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--bn")) {
            BN = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--threads")) {
            THREADS = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--log-path")) {
            LOG_PATH = arg[1];
		} else if (!strcmp(arg[0],"--save")) {
			SAVE = 1;
			arg--;
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg[0]);
            return 1;
        }
    }

	LJ_MIN = -4*LJ(R_MAX);

    if (THREADS > sysconf(_SC_NPROCESSORS_ONLN)) {
        THREADS = sysconf(_SC_NPROCESSORS_ONLN);
    }

	if (ALGO == ALGO_LISTS) {
		R = 1.2 * CUTOFF;
	} else {
		R = CUTOFF;
	}

	CSQ = CUTOFF * CUTOFF;

	if (N_PARTICLE == -1) {
		if (ALGO == ALGO_CELLS) {
			N_PARTICLE = 80 * N_CELL;
		}
		if (ALGO == ALGO_LISTS) {
			N_PARTICLE = 138 * N_CELL;
		}
	}

	printf("ALGO: %d, THREADS: %d, N_PARTICLE %d, N_CELL: %d\n", ALGO, THREADS, N_PARTICLE, N_CELL);

    return 0;
}

void init_particles(vector<particle> &particles) {
    srandom(SEED);
    for (int i = 0; i < N_PARTICLE; i++) {
        particles.push_back(particle(vec(L*frand(),L*frand(),L*frand())));
    }
}

void clipv(fvec *dst, fvec *v, float m) {
	__m256 min = _mm256_set1_ps(m);
	__m256 a = _mm256_load_ps(v);
	__m256 mask = _mm256_cmp_ps(a, min, _CMP_GT_OQ);

	a = _mm256_blendv_ps(a, min, mask);
	_mm256_store_ps(dst, a);
}

void submv(fvec *dst, fvec *va, fvec *vb) {
	__m256 a, b;
	
	a = _mm256_load_ps(va);
	b = _mm256_load_ps(vb);

	__m256 opts[3], aopts[3];
	__m256 lv = _mm256_set1_ps(L);
	opts[0] = _mm256_sub_ps(a,b);
	opts[1] = _mm256_sub_ps(opts[0],L);
	opts[1] = _mm256_add_ps(opts[0],L);

	aopts[0] = _mm256_abs_ps(opts[0]);
	aopts[1] = _mm256_abs_ps(opts[1]);
	aopts[2] = _mm256_abs_ps(opts[2]);

	__m256 m01, m12, m20;
	__m256 m0, m1, m2;
	m01 = _mm256_cmp_ps(aopts[0], aopts[1], _CMP_LT_OQ);
	m12 = _mm256_cmp_ps(aopts[1], aopts[2], _CMP_LT_OQ);
	m20 = _mm256_cmp_ps(aopts[2], aopts[0], _CMP_LT_OQ);
	
	m0 = _mm256_and_ps(m01, m12);
	m1 = _mm256_and_ps(m12, m20);
	m2 = _mm256_and_ps(m20, m01);
	
	
	__m256 res;
	res = _mm256_blendv_ps(opts[0], res, m0);
	res = _mm256_blendv_ps(opts[1], res, m1);
	res = _mm256_blendv_ps(opts[2], res, m2);

	_mm256_store_ps(dst, res);
}

void save(vector<vector<particle>> &cells, int fd) {
    if (!SAVE) 
		return;

	int ci;
    int pi;
    int np;
	particle *p;

    float buf[3];

    vector<particle> *cell;
    for (ci = 0; ci < N_CELL; ci++) {
        cell = &cells[ci];
        np = cell->size();
        for (pi = 0; pi < np; pi++) {
            p = &(*cell)[pi];
        	dprintf(fd,"%d, %f, %f, %f\n",t,p->r.x, p->r.y, p->r.z);
		}
    }
}

void save(vector<particle> &particles, int fd) {
    if (!SAVE)
		return;

	int pi;
    int np;
    particle *p;
	float buf[3];

    np = particles.size();
    
	for (pi = 0; pi < np; pi++) {
  		p = &particles[pi];
        dprintf(fd,"%d, %f, %f, %f\n", t, p->r.x, p->r.y, p->r.z);
    }
}

unsigned long rdtsc() {
  union {
    unsigned long long int64;
    struct {unsigned int lo, hi;} int32;
  } p;
  __asm__ __volatile__ (
    "rdtsc" :
    "=a" (p.int32.lo),
    "=d"(p.int32.hi)
  );
  return p.int64;
}
