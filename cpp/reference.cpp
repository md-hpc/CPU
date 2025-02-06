#include "common.h"

#include <cstdio>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <vector>

using namespace std;

char default_path[] = "particles";
char *path = default_path;

// macros for determining import regions

vector<particle> particles;


int main(int argc, char **argv) {
    char **arg;
    particle *p, p_;
    vector<particle> *hc, *nc, *cell;
    vector<particle> *outbound;

    int i, j, k, n;
    int di, dj, dk;
    int cidx, pidx, nidx;
    int ccidx[3];

    int dirfd, fd;

    float r, f, buf[3];

    vec v, v_r;

    int cur;

    
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
            exit(1);
        }
    }

    srandom(SEED);
    
    for (int i = 0; i < N_PARTICLE; i++) {
        particles.push_back(
            particle(mod_vec(L*frand(),L*frand(),L*frand()))
        );
    }

    fd = open(path,O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        perror("Could not open file");
        exit(1);
    }
    int hdr[2] = {N_PARTICLE, RESOLUTION};
    write(fd, hdr, sizeof(int) * 2);

    for (int t = 0; t < N_TIMESTEP; t++) {
        printf("Timestep %d\n",t);        
        // velocity update
        for (particle *rp = &particles[0]; rp < &particles[N_PARTICLE]; rp++) {
            for (particle *np = &particles[0]; np < &particles[N_PARTICLE]; np++) {
                if (rp->r.modr(np->r).norm() > CUTOFF)
                    continue;
                
                v = rp->r.modr(np->r);
                r = v_r.norm();
                f = lj(r);
                v *= f * DT / r;
                
                rp->v += v;
            }
        }

        // motion update
        for (particle *rp = &particles[0]; rp < &particles[N_PARTICLE]; rp++) {
            rp->r += rp->v * DT;
        }

        // save results
        if (t % RESOLUTION == 0) {
            sprintf(path,"%d",t);
            for (particle *p = &particles[0]; p < &particles[N_PARTICLE]; p++) {
                p->r.read(buf);
                write(fd,buf,3*sizeof(float));
            }
        }
    }

    close(fd);
}
