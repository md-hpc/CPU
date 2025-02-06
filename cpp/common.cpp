#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "common.h"

float SIGMA = 1;
float EPSILON = 1;
float CUTOFF = 2.5;
int UNIVERSE_SIZE = 30;
int N_PARTICLE = 2160000;
int N_TIMESTEP = 100;
float DT = 0.1;
int SEED = 0;
int RESOLUTION = 100;


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

vec &vec::operator+=(const vec& other) {
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

float vec::norm() {
    return sqrt(x*x+y*y+z*z);
}

void vec::read(float *buf) {
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
}

mod_vec::mod_vec(float x, float y, float z) : x(x), y(y), z(z) {
    normalize();
};

mod_vec mod_vec::operator+(const vec& other) {
    return mod_vec(x + other.x, y + other.y, z + other.z);
}

mod_vec &mod_vec::operator+=(const vec &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    normalize();

    return *this;
}

int mod_vec::cell() {
    int i,j,k;

    i = ((int)floorf(x/CUTOFF))%UNIVERSE_SIZE;
    j = ((int)floorf(y/CUTOFF))%UNIVERSE_SIZE;
    k = ((int)floorf(z/CUTOFF))%UNIVERSE_SIZE;
    
    return i + j * UNIVERSE_SIZE + k * UNIVERSE_SIZE * UNIVERSE_SIZE; 
}

mod_vec mod_vec::modr(mod_vec& other) {
    return mod_vec(subm(x,other.x),subm(y,other.y),subm(z,other.z));
}

void mod_vec::normalize() {
    x = fmodf(x,L);
    y = fmodf(y,L);
    z = fmodf(z,L);
}

particle::particle() : r(mod_vec(0,0,0)), v(vec(0,0,0)) {}
particle::particle(mod_vec r) : r(r), v(vec(0,0,0)) {}

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
    float c, d;
    c = b - a;
    if (c > 0) {
        d = b - L - a;
        if (c < -d) {
            return c;
        } else {
            return d;
        } 
    } else {
        d = b + L - a;
        if (-c < d) {
            return c;
        } else {
            return d;
        }
    }
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
