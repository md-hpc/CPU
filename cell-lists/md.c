#include "md.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <pthread.h>

char default_path[] = "records";

float SIGMA = 1;
float EPSILON = 1;
float CUTOFF = 2.5;
long UNIVERSE_SIZE = 3;
long N_PARTICLE = 10;
long N_TIMESTEP = 1;
float DT = 0.1;
long SEED = 0;
char *path = default_path;

long NUM_THREADS = 16;

#define N_CELL (UNIVERSE_SIZE*UNIVERSE_SIZE*UNIVERSE_SIZE)
#define L (CUTOFF * UNIVERSE_SIZE)

#define LJ_MIN (-4*24*EPSILON/SIGMA*(powf(7./26.,7./6.)-2*powf(7./26.,13./6.)))

#define DIR_MODE (S_IRWXU|S_IRWXG|S_IROTH)
#define FILE_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH)

#define particle_move(p2, p1) { \
    p2.r.x = p1.r.x; \
    p2.r.y = p1.r.y; \
    p2.r.z = p1.r.z; \
    p2.v.x = p1.v.x; \
    p2.v.y = p1.v.y; \
    p2.v.z = p1.v.z; \
}

cell_t *cells;
outbound_t *outbounds;

int main(int argc, char **argv) {
    char **arg;
    particle_t *p, _p;

    int dirfd, fd;
    pthread_t *threads;
    void *buf;

    long cidx, pidx;

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
        } else if (!strcmp(arg[0],"--threads")) {
            NUM_THREADS = atol(arg[1]);
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg[0]);
            exit(1);
        }
    }

    cells = (cell_t*) malloc(sizeof(cell_t) * N_CELL);
    for (cell_t *cell = cells; cell < &cells[N_CELL]; cell++) {
        cells->n = 0;
    }
    
    outbounds = (outbound_t*) malloc(sizeof(outbound_t) * N_CELL);

    threads = (pthread_t*) malloc(sizeof(pthread_t) * NUM_THREADS);

    srandom(SEED);
    
    p = &_p;
    for (int i = 0; i < N_PARTICLE; i++) {
        p->r.x = L*frand();
        p->r.y = L*frand();
        p->r.z = L*frand();
        
        p->v.x = 0.; //EPSILON*(frand() - .5);
        p->v.y = 0.; //EPSILON*(frand() - .5);
        p->v.z = 0.; //EPSILON*(frand() - .5);

        cidx = position_to_cell(&p->r);
        pidx = cells[cidx].n;
        particle_move(cells[cidx].particles[pidx], (*p));
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

    fd =  openat(dirfd, "spec", O_CREAT | O_TRUNC | O_RDWR, FILE_MODE);
    if (fd == -1) {
        perror("Could not open spec");
        exit(1);
    }
    dprintf(fd, "%f %f", DT, L);
    close(fd);

    for (long t = 0; t < N_TIMESTEP; t++) {
        printf("t=%d\n",t);
        // velocity update
        for (long i = 0; i < NUM_THREADS; i++) {
            pthread_create(&threads[i], NULL, velocity_update, (void*) i);
        }
        for (int i = 0; i < NUM_THREADS; i++) {
            pthread_join(threads[i], &buf);
        }

        // position update
        for (long i = 0; i < NUM_THREADS; i++) {
            pthread_create(&threads[i], NULL, position_update, (void*) i);
        }
        for (int i = 0; i < NUM_THREADS; i++) {
            pthread_join(threads[i], &buf);
        }
        
        // migration
        for (long i = 0; i < NUM_THREADS; i++) {
            pthread_create(&threads[i], NULL, particle_migration, (void*) i);
        }
        for (int i = 0; i < NUM_THREADS; i++) {
            pthread_join(threads[i], &buf);
        }
        

        // save results
        sprintf(path,"%d",t);
        fd = openat(dirfd,path,O_RDWR | O_CREAT | O_TRUNC, FILE_MODE);
        for (cell_t *cell = cells; cell < &cells[N_CELL]; cell++) {
            for (p = cell->particles; p < &cell->particles[cell->n]; p++) {
                dprintf(fd, "%f %f %f\n", p->r.x, p->r.y, p->r.z);
            }
        }
        close(fd);
    }
}

void *velocity_update(void *ptr) {
    long tid = (long) ptr;
    cell_t *hc, *nc;
    long ccidx[3];
    long hidx, nidx;
    long i, j, k;
    long di, dj, dk;
    float r, f;
    vector_t vec;

    for (hidx = tid; hidx < N_CELL; hidx += NUM_THREADS) {
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
                            if (reference == neighbor) {
                                continue;
                            }

                            modr(&vec,&reference->r,&neighbor->r);
                            r = norm(&vec);
                            
                            if (r > CUTOFF) {
                               continue;
                            }

                            if (hc != nc || vec.x < 0) {
                                continue;
                            }

                            f = lj(r);
                            scalar_mul(&vec, DT*f/r);
                            vec_add(&reference->v, &vec);
                            scalar_mul(&vec, -1.);
                            vec_add(&neighbor->v, &vec);
                        }
                    }
                }
            }
        }
    }
}

void *position_update(void *ptr) {
    long tid = (long) ptr;
    cell_t *hc;
    outbound_t *outbound;
    long cur;
    long cidx;
    particle_t *p;

    for (long hidx = tid; hidx < N_CELL; hidx+= NUM_THREADS) {
        hc = &cells[hidx];
        outbound = &outbounds[hidx];
        outbound->n = 0;
        cur = 0;
        for (particle_t *p = hc->particles; p < &hc->particles[hc->n]; p++) {
            p->r.x = fmod(p->r.x + p->v.x * DT, L);
            p->r.y = fmod(p->r.y + p->v.y * DT, L);
            p->r.z = fmod(p->r.z + p->v.z * DT, L);
            
            p->r.x = p->r.x < 0 ? p->r.x + L : p->r.x;
            p->r.y = p->r.y < 0 ? p->r.y + L : p->r.y;
            p->r.z = p->r.z < 0 ? p->r.z + L : p->r.z;

            cidx = position_to_cell(&p->r);

            if (cidx == hidx) {
                particle_move(hc->particles[cur], (*p));
                cur++;
            } else {
                particle_move(outbound->particles[outbound->n].p, (*p));
                outbound->particles[outbound->n].cidx = cidx;
                outbound->n++;
            }
        }
        hc->n = cur;
    }
}

void *particle_migration(void *ptr) {
    long tid = (long) ptr;
    cell_t *hc;
    long ccidx[3];
    long i, j, k;
    outbound_t *outbound;
    long nidx;

    for (long hidx = tid; hidx < N_CELL; hidx += NUM_THREADS) {
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
                            particle_move(hc->particles[hc->n], outbound_particle->p);
                            hc->n++;
                        }
                    }
                }
            }
        }
    }
}
 
static inline long linear_idx(long i, long j, long k) {
    i = i < 0 ? i + UNIVERSE_SIZE : i;
    j = j < 0 ? j + UNIVERSE_SIZE : j;
    k = k < 0 ? k + UNIVERSE_SIZE : k;
    return (i%UNIVERSE_SIZE) + (j%UNIVERSE_SIZE)*UNIVERSE_SIZE + (k%UNIVERSE_SIZE)*(UNIVERSE_SIZE*UNIVERSE_SIZE);
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

static inline float  subm(float a, float b) {
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
    float f = 4*EPSILON*(6*powf(SIGMA,6)/powf(r,7) - 12*powf(SIGMA,12)/powf(r,13));

    if (f < LJ_MIN) {
            return LJ_MIN;
    } else {
        return f;
    }
}

static inline float frand() {
    return ((float)rand()/(float)RAND_MAX);
}
