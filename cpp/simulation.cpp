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

	pthread_barrier_init(&barrier, NULL, THREADS);
	pthread_barrier_init(&parent_barrier, NULL, THREADS);
	tids = (pthread_t*) malloc(sizeof(pthread_t) * THREADS);
	specs = (worker_spec_t*) malloc(sizeof(worker_spec_t) * THREADS);

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
	create_workers();
	for (t = 0; t < TIMESTEPS; t++) {	
		pthread_barrier_wait(&parent_barrier);
		int np = 0;
		for (int i = 0; i < CELLS; i++)
			np += particles[i].size1();

		printf("Timestep %d, %d\n", t, np);
		if (SAVE && RESOLUTION % t == 0) {
			save();
		}
		pthread_barrier_wait(&parent_barrier);
	}
	join_workers();
}

void simulation::save() {
	int nc = particles.size();
	for (int ci = 0; ci < nc; ci++) {
		int np = particles[ci].size1();
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
	spec->s->do_work(spec);

    return (void*) 0;
}

void simulation::create_workers() {
	int t;
    int bsize =  CELLS / THREADS + (CELLS % THREADS != 0);
    
    for (t = 0; t < THREADS; t++) {
        specs[t].core = t;
        specs[t].start = t * bsize;
        specs[t].stop = (t + 1) * bsize < CELLS ? (t + 1) * bsize : CELLS;
       	specs[t].s = this;

        if (pthread_create(&tids[t], NULL, run_worker, &specs[t])) {
			perror("Could not start thread");
			exit(1);
		}
    }
}

void simulation::join_workers() {
    void *ret;
	for (t = 0; t < THREADS; t++) {
        pthread_join(tids[t], &ret);
        if (ret) {
            exit(1);
        }
    }
}
