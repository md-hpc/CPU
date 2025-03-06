#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include "simulation.h"
#include "common.h"

simulation::simulation(int argc, char **argv) {    
	THREADS = 1; // sysconf(_SC_NPROCESSORS_ONLN);
	SIGMA = 1;
	EPSILON = 1;
	CUTOFF = 1;
	DT = 1e-4;
	PARTICLES = -1;
	UNIVERSE_SIZE = 3;
	SEED = 0;
	SAVE = 0;
	TIMESTEPS = 5;
	RESOLUTION = 100;

	char *path = default_path;
	char *log_path = default_log;

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
            PARTICLES = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--timesteps")) {
            TIMESTEPS = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--dt")) {
            DT = atof(arg[1]);
        } else if (!strcmp(arg[0],"--seed")) {
            SEED = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--resolution")) {
            RESOLUTION = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--threads")) {
            THREADS = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--log-path")) {
            log_path = arg[1];
		} else if (!strcmp(arg[0],"--save")) {
			SAVE = 1;
			arg--;
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg[0]);
        }
    }

	LJ_MIN = -4*LJ(R_MAX);
	L = CUTOFF * UNIVERSE_SIZE;
	CELLS = UNIVERSE_SIZE * UNIVERSE_SIZE * UNIVERSE_SIZE;
	EP4 = 4 * EPSILON;
	SPSS = 6 * powf(SIGMA,6);
	TPST = 12 * powf(SIGMA,12);

	FD = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	LOGFD = open(log_path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

    if (THREADS > sysconf(_SC_NPROCESSORS_ONLN)) {
        THREADS = sysconf(_SC_NPROCESSORS_ONLN);
    }

	if (PARTICLES == -1)
		PARTICLES = 80 * UNIVERSE_SIZE * UNIVERSE_SIZE * UNIVERSE_SIZE;

	particles.resize(CELLS);
	outbounds.resize(CELLS);
	cios.resize(CELLS);

	srandom(SEED);
	for (int i = 0; i < PARTICLES; i++) {
		vec r = vec(L*frand(), L*frand(), L*frand());	
		vec v = vec(0,0,0);
		int hci = cell(r);
		particles[hci].append(particle(r,v,hci));
	}
}

void simulation::simulate() {
	
	for (t = 0; t < TIMESTEPS; t++) {
		int np = 0;
		for (int i = 0; i < particles.size(); i++) {
			np += particles[i].size();
		}

		printf("Timestep %d, %d\n", t, np);
		// printpv(particles);
		printf("Velocity update\n");
		thread(velocity_update, CELLS);
		// printpv(particles);

		printf("Position update\n");
		thread(position_update, CELLS);
		// printpv(particles);

		printf("Cell update\n");
		thread(cell_update, CELLS);
		// printpv(particles);

		if (SAVE && RESOLUTION % t == 0) {
			save();
		}
	}
}

void simulation::save() {
	int nc = particles.size();
	for (int ci = 0; ci < nc; ci++) {
		int np = particles[ci].size();
		for (int pi = 0; pi < np; pi++) {
			vec r = particles[ci].get(pi).r;
			printf("%d %f %f %f\n", t, r.x, r.y, r.z);
		}	
	}
}

void *simulation::run_worker(void *arg) {
    worker_spec_t *spec = (worker_spec_t*) arg;

    cpu_set_t cpuset;

    CPU_ZERO(&cpuset);
    CPU_SET(spec->core, &cpuset);
    if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
        perror("could not set affinity");
        return (void*) 1;
    }
    
    for (int i = spec->start; i < spec->stop; i++) {
        spec->worker(spec->s, i);
    }
    
    return (void*) 0;
}

void simulation::thread(worker_t worker, int n) {
    pthread_t tids[THREADS];
    worker_spec_t jobs[THREADS];
    void *ret;
    int t;
    int bsize = (n + THREADS - 1) / THREADS;
    
    for (t = 0; t < THREADS; t++) {
        jobs[t].core = t;
        jobs[t].start = t * bsize;
        jobs[t].stop = (t + 1) * bsize < n ? (t + 1) * bsize : n;
        jobs[t].worker = worker;
       	jobs[t].s = this;

        if (pthread_create(&tids[t], NULL, run_worker, &jobs[t])) {
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

void simulation::position_update(simulation *s, int hci) {
	s->position_update_worker(hci);
}

void simulation::velocity_update(simulation *s, int hci) {
	s->velocity_update_worker(hci);
}

void simulation::cell_update(simulation *s, int hci) {
	s->cell_update_worker(hci);
}
