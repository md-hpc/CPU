f8 simulation::apbcfv(f8 vx) {
	__m256 x, mlt, mgt, z, l, xu, xl;
	const float l = L;


	x = (__m256) vx;
	z = _mm256_set1_ps(0);
	l = _mm256_set1_ps(L);

	xu = (__m256) (vx + l);
	xl = (__m256) (vx - l);
	
	mlt = _mm256_cmp_ps(x,z,_CMP_LT_OQ);
	mgt = _mm256_cmp_ps(x,l,_CMP_GT_OQ);
	x = _mm256_blendv_ps(xu, x, mlt);
	x = _mm256_blendv_ps(xl, x, mgt);

	return (f8) x;
}

i8 simulation::cellv(vec8 r) {
	i8 i, j, k;
	const int u = UNIVERSE_SIZE;

	i = (i8) (r.x / CUTOFF);
	j = (i8) (r.y / CUTOFF);
	k = (i8) (r.z / CUTOFF);
	return i + j * u + k * u * u;
}

int simulation::cell(vec v) {
	const int u = UNIVERSE_SIZE;
	const int c = CUTOFF;
	int i, j, k;

	i = (int) (v.x / c);
	j = (int) (v.y / c);
	k = (int) (v.z / c);
	return i + j * u + k * u * u;
}

vec8 simulation::lj(const vec8 &rp, const vec8 &np) {
	const float ep4 = EP4;
	const float spss = SPSS;
	const float tpst = TPST;
	const float one = 1;
	const float dt = DT;
	
	f8 r = submv(rp, np);
	r = sqrtv(r.x * r.x + r.y * r.y + r.z * r.z); 

	r = one / r;
	
	f8 r8 = one, r14 = one;
	for (int i = 0; i < 8; i++)
		r8 *= r;
	for (int i = 0; i < 14; i++)
		r14 *= r;
	
	return rp * (dt * clipv(ep4 * (spss * r8 + tpst * r14), LJ_MIN));
}

vec8 simulation::submv(const vec8 &va, const vec8 &vb) {
	return vec8(
		submv(va.x,vb.x),
		submv(va.y,vb.y),
		submv(va.z,vb.z)
	);
}	

f8 simulation::submv(f8 va, f8 vb) {
	__m256 a, b;
	
	a = _mm256_load_ps(va);
	b = _mm256_load_ps(vb);

	__m256 opts[3], aopts[3];
	__m256 lv = _mm256_set1_ps(L);
	opts[0] = _mm256_sub_ps(a,b);
	opts[1] = _mm256_sub_ps(opts[0],L);
	opts[1] = _mm256_add_ps(opts[0],L);

	aopts[0] = _mm256_abs_ps(opts[0]);
	aopts[1] = _mm256_abs_ps(opts[1]);
	aopts[2] = _mm256_abs_ps(opts[2]);

	__m256 m01, m12, m20;
	__m256 m0, m1, m2;
	m01 = _mm256_cmp_ps(aopts[0], aopts[1], _CMP_LT_OQ);
	m12 = _mm256_cmp_ps(aopts[1], aopts[2], _CMP_LT_OQ);
	m20 = _mm256_cmp_ps(aopts[2], aopts[0], _CMP_LT_OQ);
	
	m0 = _mm256_and_ps(m01, m12);
	m1 = _mm256_and_ps(m12, m20);
	m2 = _mm256_and_ps(m20, m01);
	
	
	__m256 res;
	res = _mm256_blendv_ps(opts[0], res, m0);
	res = _mm256_blendv_ps(opts[1], res, m1);
	res = _mm256_blendv_ps(opts[2], res, m2);

	_mm256_store_ps(dst, res);
}

simulation::simulation(int argc, char **argv) {    
	THREADS = sysconf(_SC_NPROCESSORS_ONLN);
	SIGMA = 1;
	EPSILON = 1;
	CUTOFF = 2.5;
	DT = 1e-4;
	N_PARTICLE = -1;
	UNIVERSE_SIZE = 5;
	SEED = 0;
	SAVE = 0;
	TIMESTEPS = 20;
	RESOLUTION = 100;

	PATH = default_path;
	LOG_PATH = default_log;

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
            LOG_PATH = arg[1];
		} else if (!strcmp(arg[0],"--save")) {
			SAVE = 1;
			arg--;
        } else {
            dprintf(2,"Unrecognized option: %s\n", arg[0]);
            return 1;
        }
    }

	LJ_MIN = -4*LJ(R_MAX);
	L = CUTOFF * UNIVERSE_SIZE;
	CELLS = UNIVERSE_SIZE * UNIVERSE_SIZE * UNIVERSE_SIZE;
	EP4 = 4 * EPSILON;
	SPSS = 6 * powf(SIGMA,6);
	TPST = 12 * powf(SIGMA,12);

	FD = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	LOGFD = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

    if (THREADS > sysconf(_SC_NPROCESSORS_ONLN)) {
        THREADS = sysconf(_SC_NPROCESSORS_ONLN);
    }

	if (PARTICLES == -1)
		PARTICLES = 80 * UNIVERSE_SIZE * UNIVERSE_SIZE * UNIVERSE_SIZE;

	srandom(SEED);

	for (int i = 0; i < PARTICLES; i++) {
		vec r = vec(L*frand(), L*frand(), L*frand());	
		vec v = vec(0,0,0);
		int hci = cell(r);
		cells[hci].append(particle(r,v,hci));
	}
}

void simulation::simulate() {
	for (t = 0; t < TIMESTEPS; t++) {
		printf("Timestep %d\n", t);
		thread(velocity_update_worker, CELLS);
		thread(position_update_worker, CELLS);
		thread(cells_update_worker, CELLS);

		if (SAVE && RESOLUTION % t == 0) {
			save();
		}
	}
}

void simulation::save() {
	int nc = particles.size();
	for (int ci = 0; ci < nc; ci++) {
		int np = particles[i].size();
		for (int pi = 0; pi < np; pi++) {
			vec r = particles[i].get(pi).r;
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
        spec->worker(i);
    }
    
    return (void*) 0;
}

void simulation::thread(void (*worker)(int), int n) {
    pthread_t tids[THREADS];
    job_t jobs[THREADS];
    void *ret;
    int t;
    int bsize = (n + THREADS - 1) / THREADS;
    
    for (t = 0; t < THREADS; t++) {
        jobs[t].core = t;
        jobs[t].start = t * bsize;
        jobs[t].stop = (t + 1) * bsize < n ? (t + 1) * bsize : n;
        jobs[t].worker = worker;
        
        if (pthread_create(&tids[t], NULL, run_kernel, &jobs[t])) {
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

voxel simuation::voxelof(int idx) {
	int i,j,k;
	const int u = UNIVERSE_SIZE;

	k = idx / (u * u);
    idx -= res[2] * u * u;
    j = idx / u;
    idx -= res[1] * u;
    i = idx;
	
	return voxel(i,j,k);
}

int simulation::cell(int i, int j, int k) {
	const int u = UNIVESE_SIZE;
	i = i < 0 ? i+u : i;
	j = j < 0 ? j+u : j;
	k = k < 0 ? k+u : k;

	return i + u * j + u * u * k;
}

voxel::voxel(int i, int j, int k) : i(i), j(j), k(k) {}
