#include <stdio.h>
#include <math.h>

#include "simulation.h"
#include "avx.h"
#include "particle8.h"
#include "common.h"

void simulation::do_work(worker_spec_t *spec) {
	int start = spec->start;
	int stop = spec->stop;

	for (int t = 0; t < TIMESTEPS; t++) {
		pthread_barrier_wait(&parent_barrier);
		// let the parent do its thing
		pthread_barrier_wait(&parent_barrier);
	
		for (int i = start; i < stop; i++) {
			velocity_update_worker(i,spec);
		}
		pthread_barrier_wait(&barrier);

		for (int i = start; i < stop; i++) {
			position_update_worker(i,spec);
		}
		pthread_barrier_wait(&barrier);

		for (int i = start; i < stop; i++) {
			cell_update_worker(i,spec);
		}
	}
}

void simulation::velocity_update_worker(int hci, worker_spec_t* spec) {
	voxel hcv = voxelof(hci);
	int nr = particles[hci].size8();
	
	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				if (di < 0 || di == 0 && dj < 0 || di == 0 && dj == 0 && dk < 0)
					continue;

				int nci = cell(hcv.i + di, hcv.j + dj, hcv.k + dk);				
				int nn = particles[nci].size8();
				
				for (int ri = 0; ri < nr; ri++) {
					particle8 rp = particles[hci][ri];

					for (int p = 0; p < VSIZE; p++) {
						for (int ni = 0; ni < nn; ni++) {
							// nci >= hci from conditional before
							//
							// ni > ri (&& nci == hci) should be omitted to leverag n3l
							// ni == ri (&& nci == hci) requires special handling (no n3l)
							if (nci == hci && ni >= ri) 
								continue;

							particle8 np = particles[nci][ni];

							vec8 v = lj(rp.r, np.r);
							rp.v += v;

							v *= -1;
							particles[nci][ni].v += v;
						}
				
						rp.r.permutev();
						rp.v.permutev();
					}
					particles[hci][ri].v = rp.v;
				}
			}
		}
	}

	// handle ri == ni && hci == nci for this cell. Do not apply n3l
	for (int i = 0; i < nr; i++) {
		vec8 rp, np, rv, v;
		rp = particles[hci][i].r;
		rv = particles[hci][i].v;
		np = particles[hci][i].r;

		for (int p = 0; p < 7; p++) {
			np.permutev();
			vec8 v = lj(rp, np);
			rv += v;
		}
		particles[hci][i].v = rv;
	}
}

void simulation::position_update_worker(int hci, worker_spec_t *spec) {
	// hci === home cell index

	int np = particles[hci].size8();
	int cur = 0;
	const int hciv = hci;

	// buffer for particles that have left this cell
	outbounds[hci].resize(0);

	// consolidation buffer
	p8buf buf;

	for (int pi = 0; pi < np; pi++) {
		particle8 p = particles[hci][pi];
		p.r += (p.v * DT);
		apbcfv(p.r); // apply periodic boundary condition (floating-point, vector)

		ipack cells = { .v=cellv(p.r) };

		if (alleq(cells.v, hciv)) {
			// if no particles have left this cell, perform aligned move
			particles[hci][cur++] = p;
		} else {
			// we must pick particle-by-particle which need to be moved to the outbound buffer
			for (int i = 0; i < VSIZE; i++) {
				if (cells.d[i] == hci) {
					// append to consolidation buffer
					if (buf.append(p.get(i))) {
						// if buffer is full, perform aligned store to the cell list
						particles[hci][cur++] = buf.get();
					}
				} else if (cells.d[i] > -1) {
					// append to outbound buffer
					// 
					// if this is a nonexistent particle (NAN) its cell will be a large negative number, so
					// we filter out those with the conditional clause
					
					// TODO: if op.cell is in this thread's group, then just
					// append it to that cell's list
					particle op = p.get(i);
					op.cell = cells.d[i];
					outbounds[hci].push_back(op);
				}
			}
		}
	}
	

	int sz = cur * VSIZE;
	if (buf.i > 0) {
		// flush remaining buffer to the cell list
		sz += buf.i;
		particles[hci][cur++] = buf.get();
	}
	particles[hci].resize(sz);
}

void simulation::cell_update_worker(int hci, worker_spec_t *spec) {
	// hci == home cell index
	//
	// check all neighbor outbound buffers to see if they belong to this cell
	voxel hcv = voxelof(hci);
	
	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				int oci = cell(hcv.i + di, hcv.j + dj, hcv.k + dk);
				int no = outbounds[oci].size();	
				for (int oi = 0; oi < no; oi++) {
					particle op = outbounds[oci][oi];
					if (op.cell == hci) {
						particles[hci].append(op);
					}
				}
			}
		}
	}
}
