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
char default_log[] = "validate/reference.out";
char *path = default_path;
char *LOG_PATH = default_log;

// macros for determining import regions

vector<particle> particles;

int t;

int main(int argc, char **argv) {
    char **arg;

	ALGO = ALGO_CELLS;
    parse_cli(argc,argv);
    
    init_particles(particles);

	int log_fd = open(LOG_PATH, O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR);
	if (log_fd == -1) {
		perror("Could not open file");
		exit(1);
	}

    for (t = 0; t < N_TIMESTEP; t++) {
        printf("Timestep %d\n",t);
        // velocity update
        particle *pr, *pn;
        int ri, ni;
        for (ri = 0; ri < N_PARTICLE; ri++) {
            pr = &particles[ri];
            for (ni = 0; ni < N_PARTICLE; ni++) {
                pn = &particles[ni];

				if (pr->interact(pn)) {
					#ifdef DEBUG 
					dprintf(log_fd,"%d %d %d\n",t,pr->id,pn->id);
					#endif	
				}
				
				#ifdef DEBUG
                if (pr->id == BR && pn->id == BN)
                    printf("%d %d\n", BR, BN);    
				if (pr->id == BR)
					printf("%d\n", BR);
				#endif
            }
        }

        // position update
		for (int pi = 0; pi < N_PARTICLE; pi++) {
            particles[pi].update_position();
        }

        // save result 
        if (t % RESOLUTION == 0) {
            save(particles, log_fd);
        }
    }
}
