// cell lists using hemisphere method

#include "common.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <vector>

using namespace std;

char default_path[] = "hemisphere-particles";

char *path = default_path;

vector<vector<particle*>> neighbors;
vector<vector<particle>> cells;
vector<vector<particle>> outbounds;


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

    parse_cli(argc, argv);

    cells.resize(N_CELL);
    outbounds.resize(N_CELL); 

    srandom(SEED); 
    for (int i = 0; i < N_PARTICLE; i++) {
        p_ = particle(mod_vec(L*frand(),L*frand(),L*frand()));
        cells[p_.cell()].push_back(p_);
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
        if (t % NEIGHBOR_REFRESH_RATE == 0) {
             
            for (int hidx = 0; hidx < N_CELL; hidx++) {
                hc = &cells[hidx];

                cubic_idx(ccidx, hidx);
                i = ccidx[0];
                j = ccidx[1];
                k = ccidx[2];

                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int dk = -1; dk <= 1; dk++) {
                            if (di < 0 || di == 0 && dj < 0 || di < 0 && dj < 0 && dk < 0) {
                                continue;
                            }
                            nidx = linear_idx(i+di, j+dj, k+dk);
                            nc = &cells[nidx];
                            
                            particle *pr, *pn;
                            int nr, nn; 
                            for (pr = &(*hc)[0], nr = hc->size(); pr < &(*hc)[nr]; pr++) {
                                for (pn = &(*nc)[0], nn = nc->size(); pn < &(*nc)[nn]; pn++) {
                                    if (pr->r.x > pn->r.x)
                                        continue;

                                    v_r = pr->r.modr(pn->r);            
                                    r = v_r.norm();

                                    if (r > CUTOFF || r == 0) {
                                       continue;
                                    }
                                    f = lj(r);
                                    v_r *= DT*f/r;
                                    
                                    pr->v += v_r;
                                    v_r *= -1.;
                                    pn->v += v_r;
                                }
                            }
                        }
                    }
                }
            }
        
        // position update
        for (int i = 0; i < outbounds.size(); i++) {
            outbounds[i].resize(0);
        }
        
        int hidx, nc, np;
        vector<particle> *hc;
        particle *p;
        for (hidx = 0, nc = cells.size(), hc = &cells[hidx]; hidx < nc; hc = &cells[++hidx]) {
            cur = 0;
            
            for (p = &(*hc)[0], np = hc->size(); p < &(*hc)[np]; p++) {
                p->r += p->v * DT;
                
                cidx = p->cell();
                if (cidx == hidx) {
                    (*hc)[cur++] = *p;
                } else {
                    outbounds[hidx].push_back(*p);
                }
            }
            hc->resize(cur);
        }
        
        // particle migration
        for (int hidx = 0; hidx < N_CELL; hidx++) {
            hc = &cells[hidx];
            cubic_idx(ccidx, hidx);
            i = ccidx[0];
            j = ccidx[1];
            k = ccidx[2];

            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    for (int dk = -1; dk <= 1; dk++) {
                        nidx = linear_idx(i+di, j+dj, k+dk);
                        outbound = &outbounds[nidx];
                        
                        particle *p;
                        int np;
                        for (p = &(*outbound)[0], np = outbound->size(); p < &(*outbound)[np]; p++) {
                            if (p->cell() == hidx) {
                                hc->push_back(*p);
                            }
                        }
                    }
                }
            }
        }

        // save results

        if (t % RESOLUTION == 0) {
            sprintf(path,"%d",t);
            for (vector<vector<particle>>::iterator cell = cells.begin(); cell != cells.end(); ++cell) {
                for (vector<particle>::iterator p = cell->begin(); p != cell->end(); ++p) {
                    p->r.read(buf);
                    write(fd,buf,3*sizeof(float));
                }
            }
        }
    }

    close(fd);
}
