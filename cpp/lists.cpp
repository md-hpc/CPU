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

char default_path[] = "particles";
char default_log[] = "validate/lists/interactions";
char *LOG_PATH = default_log;
char *path = default_path;


void make_neighbor_lists(int,int);
void velocity_update(int,int);
void position_update(int,int);


vector<vector<vector<particle*>>> neighbors;
vector<vector<particle>> cells;

vector<particle> particles;

int t;

int *fds;

int main(int argc, char **argv) {
    particle *p;
	int fd;

	ALGO = ALGO_LISTS;
	parse_cli(argc, argv);

    cells.resize(N_CELL);
    neighbors.resize(N_CELL);
    
    init_particles(particles);
    for (int i = 0; i < N_PARTICLE; i++) {
        p = &particles[i];
        cells[p->cell].push_back(*p);
    }

	fd = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

    #ifdef DEBUG
	fds = (int*) malloc(sizeof(int) * THREADS);
    for (int tid = 0; tid < THREADS; tid++) {
		char path[16];
		sprintf(path,"%s%d",LOG_PATH,tid);
		fds[tid] = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	}
	#endif

    for (t = 0; t < N_TIMESTEP; t++) {
        printf("Timestep %d\n",t);      
        
        // velocity update
        if (t % NEIGHBOR_REFRESH_RATE == 0) {
            // collect particles into one big list
            int cidx;

            for (cidx = 0; cidx < N_CELL; cidx++) {
                vector<particle> *cell = &cells[cidx];
                int np = cell->size();
                for (int pidx = 0; pidx < np; pidx++) {
                    particles.push_back((*cell)[pidx]);
                }
                cell->resize(0);
            }

            // bin particles into cells
            for (int pidx = 0; pidx < N_PARTICLE; pidx++) {
                cidx = particles[pidx].cell;
                cells[cidx].push_back(particles[pidx]);
            }
            
            // resize 2D neighbor lists
            for (int cidx = 0; cidx < N_CELL; cidx++) {
                neighbors[cidx].resize(cells[cidx].size());
            }

            // construct neighbor lists
            thread(make_neighbor_lists, N_CELL);
        }

        // velocity update
        thread(velocity_update, N_CELL);

        // position update       
        thread(position_update, N_CELL);

        if (t % RESOLUTION == 0) {
        	save(cells,fd);
		}
    }

    close(fd);
}

void make_neighbor_lists(int hci, int tid) {
    int ccidx[3];
    float r;
    int i,j,k;
	vec v;

    vector<particle> *hc = &cells[hci];
    cubic_idx(ccidx, hci);
    i = ccidx[0];
    j = ccidx[1];
    k = ccidx[2];
    
	float csq = R * R;

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
                int nci = linear_idx(i+di, j+dj, k+dk);
                vector<particle> *nc = &cells[nci];
                
                particle *pr, *pn;
                int ri, ni;
                int nr, nn; 
                for (ri = 0, nr = hc->size(); ri < nr; ri++) {
                    pr = &cells[hci][ri];
                    vector<particle*> *neighbor_list = &neighbors[hci][ri];

                    #ifdef DEBUG
                    if (pr->id == BR) {
                        printf("%d\n",BR);
                    }
                    #endif

                    for (ni = 0, nn = nc->size(); ni < nn; ni++) {
                        pn = &cells[nci][ni];

                       	if (pn == pr)
							continue;
                        
						#ifdef DEBUG
                        if (pn->id == BN && pr->id == BR) {
                            printf("%d %d\n",BR, BN);
                        }
                        #endif

                        v = pn->r % pr->r;
						if (v.x < 0)
							continue;

						r = v.normsq();
                        if (r < csq) {
                            neighbor_list->push_back(pn);
                        }
                    }
                }
            }
        }
    }
}

void velocity_update(int ci, int tid) {
    vec v;
    float r,f;
	float csq = R * R;

    vector<particle> *cell = &cells[ci];
    vector<vector<particle*>> *cell_neighbors = &neighbors[ci];
    int nr = cell_neighbors->size();
    for (int ri = 0; ri < nr; ri++) { 
        particle *pr = &cells[ci][ri];

        vector<particle*> *neighbor_list = &(*cell_neighbors)[ri];
        int nn = neighbor_list->size();
        for (int ni = 0; ni < nn; ni++) {
            particle *pn = (*neighbor_list)[ni];
                                
			#ifdef DEBUG
            if (pr->id == BR && pn->id == BN) {
                printf("%d %d\n", BR, BN);
            }
			#endif

            if (pr->interact(pn)) {
				#ifdef DEBUG
				dprintf(fds[tid],"%d %d %d\n", t, pr->id, pn->id);
				#endif
			}
        }
    }
}

void position_update(int ci, int tid) {
    vector<particle> *cell = &cells[ci];
    int np = cell->size();
    particle *p;

    for (int pi = 0; pi < np; pi++) {
        particles[pi].update_position();
    }
}
