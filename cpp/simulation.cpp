fvec simulation::apbcfv(fvec vx) {
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

	return (fvec) x;
}

ivec simulation::cellv(fvec x, fvec y, fvec z) {
	ivec i, j, k;
	const int u = UNIVERSE_SIZE;

	i = (ivec) (x / CUTOFF);
	j = (ivec) (j / CUTOFF);
	k = (ivec) (k / CUTOFF);
	return i + j * u + k * u * u;
}

fvec simulation::ljv(fvec r) {
	const float ep4 = 4 * EPSILON;
	const float spss = 6 * powf(SIGMA,6);
	const float tpst = 12 * powf(SIGMA,12);
	const float one = 1;

	fvec r = sqrtv(r); 

	r = one / r;
	
	fvec r8 = 1, r14 = 1;
	for (int i = 0; i < 8; i++)
		r8 *= r;
	for (int i = 0; i < 14; i++)
		r14 *= r;
	
	return clipv(ep4 * (spss * r8 + tpst * r14), LJ_MIN);
}


fvec simulation::submv(fvec va, fvec vb) {
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

int simulation::parse_cli(int argc, char **argv) {    
    if (ALGO == ALGO_NONE) {
		printf("Must set ALGO before calling parse_cli");
		return 1;
	}

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
        } else if (!strcmp(arg[0],"--br")) {
            BR = atoi(arg[1]);
        } else if (!strcmp(arg[0],"--bn")) {
            BN = atoi(arg[1]);
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

    if (THREADS > sysconf(_SC_NPROCESSORS_ONLN)) {
        THREADS = sysconf(_SC_NPROCESSORS_ONLN);
    }

	if (ALGO == ALGO_LISTS) {
		R = 1.2 * CUTOFF;
	} else {
		R = CUTOFF;
	}

	CSQ = CUTOFF * CUTOFF;

	if (N_PARTICLE == -1) {
		if (ALGO == ALGO_CELLS) {
			N_PARTICLE = 80 * N_CELL;
		}
		if (ALGO == ALGO_LISTS) {
			N_PARTICLE = 138 * N_CELL;
		}
	}

	printf("ALGO: %d, THREADS: %d, N_PARTICLE %d, N_CELL: %d\n", ALGO, THREADS, N_PARTICLE, N_CELL);

    return 0;
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
