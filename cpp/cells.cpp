// cell lists. Do not apply N3L

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
char default_log[] = "validate/cells/interactions";
char *LOG_PATH = default_log;


// holds particles that have exceeded the bounds of their cells, but have yet to be migrated to their new cell
// we'll need this for thread-safety as there cannot be multiple writers to an std::vector
vector<vector<particle>> outbounds;
vector<vector<particle>> cells;
vector<particle> particles;

void velocity_update(int,int);
void position_update(int,int);
void cell_update(int,int);

int t;

int *fds;

int main(int argc, char **argv) {
    char **arg;
    particle *p, p_;
    vector<particle> *hc, *nc, *cell;
    vector<particle> *outbound;

    int i, j, k, n;
    int di, dj, dk;
    int ci, pi, nci;
    int cci[3];

    int dirfd, fd;

    float r, f, buf[3];

    vec v, v_r;

    int cur;

	ALGO = ALGO_CELLS;
    if (parse_cli(argc, argv)) {
        exit(1);
    }

    cells.resize(N_CELL);
    outbounds.resize(N_CELL); 

    init_particles(particles);
    for (int i = 0; i < N_PARTICLE; i++) {
        p = &particles[i];
        cells[p->cell].push_back(*p);
    }
    particles.resize(0);


	fd = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

	#ifdef DEBUG
	fds = (int*) malloc(sizeof(int) * THREADS);
	for (int tid = 0; tid < THREADS; tid++) {
		char path[16]; 
		sprintf(path,"%s%d", LOG_PATH, tid); 
    	fds[tid] = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		if (fds[tid] == -1) {
        	perror("Could not open file");
        	exit(1);
    	}
	}
	#endif

    for (t = 0; t < N_TIMESTEP; t++) {
        printf("Timestep %d\n",t);
        
        thread(velocity_update, N_CELL);
 
        thread(position_update, N_CELL);

        thread(cell_update, N_CELL);
        
        // save results
        if (t % RESOLUTION == 0) {
            save(cells, fd);
        }
    }

    close(fd);
}

void velocity_update(int hci, int tid) {
    int i,j,k;
    vector<particle> *hc, *nc;
    int cci[3];
    int nci;
    float r, f;
    vec v;

    hc = &cells[hci];

    cubic_idx(cci, hci);
    i = cci[0];
    j = cci[1];
    k = cci[2];

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {

				nci = linear_idx(i+di, j+dj, k+dk);
                nc = &cells[nci];
                
                particle *pr, *pn;
                int ri, ni;
                int nr = hc->size(), nn = nc->size(); 
                for (ri = 0; ri < nr; ri++) {
                    pr = &(*hc)[ri];
                    for (ni = 0; ni < nn; ni++) {
                        pn = &(*nc)[ni];
						
						if (pr->interact(pn)) {
							#ifdef DEBUG
                        	dprintf(fds[tid], "%d %d %d\n", t, pr->id, pn->id);
							#endif
						}

						
						#ifdef DEBUG
                        if (pr->id == BR && pn->id == BN) {
                            printf("%d %d\n");
                        }
						#endif
                    }
                }
            }
        }
    }
}

void position_update(int ci, int tid) {
    int cur, np, hci, pi;
    particle *p;
    vector<particle> *cell;

    outbounds[ci].resize(0);
    cell = &cells[ci];
    cur = 0;
    np = cell->size();
    for (pi = 0; pi < np; pi++) {
        p = &(*cell)[pi];
		p->update_position();

        hci = p->cell;
        if (hci == ci) {
            (*cell)[cur++] = *p;
        } else {
            outbounds[ci].push_back(*p);
        }
    }
    cell->resize(cur);
}

void cell_update(int hci, int tid) {
    vector<particle> *hc, *outbound;
    int i, j, k;
    int nci;
    int cci[3];

    hc = &cells[hci];
    cubic_idx(cci, hci);
    i = cci[0];
    j = cci[1];
    k = cci[2];

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
                nci = linear_idx(i+di, j+dj, k+dk);
                outbound = &outbounds[nci];
                
                int pi;
                particle *p;
                int np = outbound->size();
                for (pi = 0; pi < np; pi++) {
                    p = &(*outbound)[pi];
                    if (p->cell == hci) {
                        hc->push_back(*p);
                    }
                }
            }
        }
    }
}
