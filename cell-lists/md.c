#include "md.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>

char default_path[] = "particles";

float SIGMA = 1;
float EPSILON = 1;
float CUTOFF = 2.5;
long UNIVERSE_SIZE = 3;
long N_PARTICLE = 10;
long N_TIMESTEP = 5;
float DT = 0.1;
long SEED = 0;
char *path = default_path;

#define N_CELL (UNIVERSE_SIZE*UNIVERSE_SIZE*UNIVERSE_SIZE)
#define L (CUTOFF * UNIVERSE_SIZE)

#define LJ_MIN (-4*24*EPSILON/SIGMA*(powf(7./26.,7./6.)-2*powf(7./26.,13./6.)))

#define DIR_MODE (S_IRWXU|S_IRWXG|S_IROTH)
#define FILE_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH)

cell_t *cells;

int main(int argc, char **argv) {
    char **arg;
    particle_t *p, _p;
    cell_t *hc, *nc, *cell;

    long i, j, k;
    long di, dj, dk;
    long cidx, pidx, nidx;
    long ccidx[3];

    int dirfd, fd;

    vector_t vec;
    float r, f;

    long cur;
    outbound_t *outbound, *outbounds;

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
            N_PARTICLE = atol(arg[1]);
        } else if (!strcmp(arg[0],"--timesteps")) {
            N_TIMESTEP = atol(arg[1]);
        } else if (!strcmp(arg[0],"--dt")) {
            DT = atof(arg[1]);
        } else if (!strcmp(arg[0],"--seed")) {
            SEED = atol(arg[1]);
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg);
            exit(1);
        }
    }

    cells = (cell_t*) malloc(sizeof(cell_t) * N_CELL);
    for (cell_t *cell = cells; cell < &cells[N_CELL]; cell++) {
        cells->n = 0;
    }
    
    outbounds = (outbound_t*) malloc(sizeof(outbound_t) * N_CELL);

    srandom(SEED);
    
    p = &_p;
    for (int i = 0; i < N_PARTICLE; i++) {
        p->r.x = L*frand();
        p->r.y = L*frand();
        p->r.z = L*frand();
        
        p->v.x = EPSILON*(frand() - .5);
        p->v.y = EPSILON*(frand() - .5);
        p->v.z = EPSILON*(frand() - .5);

        cidx = position_to_cell(&p->r);
        pidx = cells[cidx].n;
        memcpy(&cells[cidx].particles[pidx], p, sizeof(particle_t));
        cells[cidx].n++;
    }


    if (access(path, F_OK) && mkdir(path, DIR_MODE)) {
        perror("Could not create directory");
        exit(1);
    }
    dirfd = open(path, O_DIRECTORY | O_RDONLY);
    if (dirfd == -1) {
        perror("Could not open directory");
        exit(1);
    }


    for (long t = 0; t < N_TIMESTEP; t++) {
        // velocity update
        for (long hidx = 0; hidx < N_CELL; hidx++) {
            hc = &cells[hidx];

            cubic_idx(ccidx, hidx);
            i = ccidx[0];
            j = ccidx[1];
            k = ccidx[2];

            for (long di = -1; di <= 1; di++) {
                for (long dj = -1; dj <= 1; dj++) {
                    for (long dk = -1; dk <= 1; dk++) {
                        if (di < 0 || di == 0 && dj < 0 || di < 0 && dj < 0 && dk < 0) {
                            continue;
                        }
                        nidx = linear_idx(i+di, j+dj, k+dk);
                        nc = &cells[nidx];

                        for (particle_t *reference = hc->particles; reference < &hc->particles[hc->n]; reference++) {
                            for (particle_t *neighbor = nc->particles; neighbor < &nc->particles[nc->n]; neighbor++) {
                                modr(&vec,&reference->r,&neighbor->r);
                                r = norm(&vec);
                                
                                if (r > CUTOFF || r == 0) {
                                   continue;
                                }
                                f = lj(r);
                                scalar_mul(&vec, DT*f/r);
                                vec_add(&reference->v, &vec);
                                scalar_mul(&vec, -1.);
                                neighbor->v.x += vec.x;
                                neighbor->v.y += vec.y;
                                neighbor->v.z += vec.z;
                            }
                        }
                    }
                }
            }
        }
        
        // position update
        for (long hidx = 0; hidx < N_CELL; hidx++) {
            hc = &cells[hidx];
            outbound = &outbounds[hidx];
            outbound->n = 0;
            cur = 0;

            for (particle_t *p = hc->particles; p < &hc->particles[hc->n]; p++) {
                p->r.x = fmod(p->r.x + p->v.x * DT, L);
                p->r.y = fmod(p->r.y + p->v.y * DT, L);
                p->r.z = fmod(p->r.z + p->v.z * DT, L);
                
                cidx = position_to_cell(&p->r);
                if (cidx == hidx) {
                    memmove(&hc->particles[cur], p, sizeof(particle_t));
                    cur += 1;
                } else {
                    memcpy(&outbound->particles[outbound->n].p, p, sizeof(particle_t));
                    outbound->particles[outbound->n].cidx = cidx;
                    outbound->n++;
                }
            }
            hc->n = cur;
        }
        
        // particle migration
        for (long hidx = 0; hidx < N_CELL; hidx++) {
            hc = &cells[hidx];
            cubic_idx(ccidx, hidx);
            i = ccidx[0];
            j = ccidx[1];
            k = ccidx[2];

            for (long di = -1; di <= 1; di++) {
                for (long dj = -1; dj <= 1; dj++) {
                    for (long dk = -1; dk <= 1; dk++) {
                        nidx = linear_idx(i+di, j+dj, k+dk);
                        outbound = &outbounds[nidx];
                        for (outbound_particle_t *outbound_particle = outbound->particles; outbound_particle < &outbound->particles[outbound->n]; outbound_particle++) {
                            if (outbound_particle->cidx == hidx) {
                                memcpy(&hc->particles[hc->n], &outbound_particle->p, sizeof(particle_t));
                                hc->n++;
                            }
                        }
                    }
                }
            }
        }

        // save results
        sprintf(path,"%d",t);
        fd = openat(dirfd,path,O_RDWR | O_CREAT | O_TRUNC);
        for (cell = cells; cell < &cells[N_CELL]; cell++) {
            for (p = cell->particles; p < &cell->particles[cell->n]; p++) {
                write(fd,&p->r,sizeof(vector_t));
            }
        }
        close(fd);
    }
}


static inline long linear_idx(long i, long j, long k) {
    i = i < 0 ? i + UNIVERSE_SIZE : i;
    j = j < 0 ? j + UNIVERSE_SIZE : j;
    k = k < 0 ? k + UNIVERSE_SIZE : k;
    return (i%UNIVERSE_SIZE) + (j%UNIVERSE_SIZE)*UNIVERSE_SIZE + (k%UNIVERSE_SIZE)*powl(UNIVERSE_SIZE,2);
}

static inline void cubic_idx(long *res, long idx) {
    res[0] = idx % UNIVERSE_SIZE;
    res[1] = (idx/UNIVERSE_SIZE)%UNIVERSE_SIZE;
    res[2] = (idx/(UNIVERSE_SIZE*UNIVERSE_SIZE))%UNIVERSE_SIZE;
}

static inline long position_to_cell(vector_t *r) {
    long i, j, k;

    i = (long)floorf(r->x/CUTOFF)%UNIVERSE_SIZE;
    j = (long)floorf(r->y/CUTOFF)%UNIVERSE_SIZE;
    k = (long)floorf(r->z/CUTOFF)%UNIVERSE_SIZE;
    
    return linear_idx(i, j, k);
}

static inline float subm(float a, float b) {
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

// these can be vectorized with SIMD
static inline void scalar_mul(vector_t *v, float c) {
    v->x *= c;
    v->y *= c;
    v->z *= c;
}

static inline void vec_add(vector_t *a, vector_t *b) {
    a->x += b->x;
    a->y += b->y;
    a->z += b->z;
}

static inline void modr(vector_t *c, vector_t *a, vector_t *b) {
    c->x = subm(a->x, b->x);
    c->y = subm(a->y, b->y);
    c->z = subm(a->z, b->z);
}

static inline float norm(vector_t *r) {
    return sqrt(powf(r->x,2) + powf(r->y,2) + powf(r->z,2)); 
}

static inline float lj(float r) {
    float f = 4*EPSILON*(6*pow(SIGMA,6)/pow(r,7) - 12*pow(SIGMA,12)/pow(r,13));

    if (f < LJ_MIN) {
            return LJ_MIN;
    } else {
        return f;
    }
}

static inline float frand() {
    return ((float)rand()/(float)RAND_MAX);
}
