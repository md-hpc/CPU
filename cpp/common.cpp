#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "common.h"
#include <unistd.h>

using namespace std;

float SIGMA = 1;
float EPSILON = 1;
float CUTOFF = 2.5;
int UNIVERSE_SIZE = 3;
int N_PARTICLE = 2160;
int N_TIMESTEP = 10;
float DT = 0.1;
int SEED = 0;
int RESOLUTION = 100;
int NEIGHBOR_REFRESH_RATE = 10;


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
}

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

void vec::apbc() {
    x = fmodf(x,L);
    y = fmodf(y,L);
    z = fmodf(z,L);
}

int vec::cell() {
    return linear_idx(
        (int)x/CUTOFF,
        (int)y/CUTOFF,
        (int)z/CUTOFF
    );
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
    id = counter++;
}

int particle::cell() {
    return r.cell();
}

int linear_idx(int i, int j, int k) {
    i = i < 0 ? i + UNIVERSE_SIZE : i;
    j = j < 0 ? j + UNIVERSE_SIZE : j;
    k = k < 0 ? k + UNIVERSE_SIZE : k;
    return (i%UNIVERSE_SIZE) + (j%UNIVERSE_SIZE)*UNIVERSE_SIZE + (k%UNIVERSE_SIZE)*powl(UNIVERSE_SIZE,2);
}

void cubic_idx(int *res, int idx) {
    res[0] = idx % UNIVERSE_SIZE;
    res[1] = (idx/UNIVERSE_SIZE)%UNIVERSE_SIZE;
    res[2] = (idx/(UNIVERSE_SIZE*UNIVERSE_SIZE))%UNIVERSE_SIZE;
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
    float f = 4*EPSILON/r*(6*pow(SIGMA,6)/pow(r,7) - 12*pow(SIGMA,12)/pow(r,13));

    if (f < LJ_MIN) {
            return LJ_MIN;
    } else {
        return f;
    }
}

float frand() {
    return ((float)rand()/(float)RAND_MAX);
}


void thread(void (*kernel)(int), long n) {
    pthread_t tids[n];
    long i;
    void *ret;

    for (i = 0; i < n; i++) {
        pthread_create(&tids[i], NULL, (void*(*)(void*)) kernel, (void *) i);
    }

    for (i = 0; i < n; i++) {
        pthread_join(tids[i], &ret);
    }
}

int parse_cli(int argc, char **argv) {    
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
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg);
            return 1;
        }
    }

    return 0;
}

void init_particles(vector<particle> &particles) {
    srandom(SEED);
    for (int i = 0; i < N_PARTICLE; i++) {
        particles.push_back(particle(vec(L*frand(),L*frand(),L*frand())));
    }
}

void save(vector<vector<particle>> &cells, int fd) {
    int ci;
    int pi;
    int np;

    float buf[3];

    vector<particle> *cell;
    for (ci = 0; ci < N_CELL; ci++) {
        cell = &cells[ci];
        np = cell->size();
        for (pi = 0; pi < np; pi++) {
            (*cell)[pi].r.read(buf);
            write(fd,buf,3*sizeof(float));
        }
    }
}

void save(vector<particle> &particles, int fd) {
    int pi;
    int np;
    float buf[3];

    np = particles.size();
    for (pi = 0; pi < np; pi++) {
        particles[pi].r.read(buf);
        write(fd,buf,3*sizeof(float));
    }
}
