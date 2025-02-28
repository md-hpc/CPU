#include "simulation.h"
#include "avx.h"
#include "fvv.h"
#include "common.h"

void simulation::velocity_update_worker(int hci) {
	voxel hcv = voxelof(hci);
	int nr = particles[hci].size();

	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				if (di < 0 || di == 0 && dj < 0 || di == 0 && dj == 0 && dk < 0)
					continue;

				int nci = cell(hcv.i + di, hcv.j + dj, hcv.k + dk);				
				int nn = particles[nci].size();
				
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

							vec8 v = ljv(rp.r, np.r);
							rp.v += v;

							v *= -1;
							particles[nci][ni].v += v;
						}
				
						rp.r = permute(rp.r);
						rp.v = permute(rp.v);
					}
					v[hci][ri].v = rp.v;
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
			np = permute(np);
			vec8 v = ljv(rp, np);
			rv += v;
		}
		particles[hci][i].v = rv;
	}
}

void simulation::position_update_worker(int hci) {

	int np = r[hci].size();
	int cur = 0;
	int ocur = 0;

	outbounds[hci].resize(0);

	for (int pi = 0; pi < np; pi++) {
		particles[hci][pi] = apbcfv(
			particles[hci][pi].r + particles[hci][pi].r * DT
		);

		ipack cells = { .v=cellv(r[hci].x[pi], r[hci].y[pi], r[hci].z[pi]) };

		// I think this could be optimized if we tried to do as many aligned vector loads as possible when
		// all of cells == hci
		for (int i = pi * VSIZE; i < (pi + 1) * VSIZE; i++) {
			if (cells.d[i%VSIZE] == hci) {
				r[hci].set(r[hci].get(i),cur);
				v[hci].set(v[hci].get(i),cur);
				cur++;
			} else {
				ros[hci].push_back(r[hci].get(i));
				vos[hci].push_back(v[hci].get(i));
				cios[hci].push_back(cells.d[i%VSIZE]);
			}
		}
	}
	v[hci].resize(cur);
	r[hci].resize(cur);
}

void simulation::cell_update_worker(int hci) {
	voxel hcv = voxelof(hci);
	
	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				int nci = cell(hcv.i + di, hcv.j + dj, hcv.k + dk);
				int no = cios[nci].size();
				for (int oi = 0; oi < no; oi++) {
					if (cios[nci][oi] == hci) {
						r[hci].append(ros[nci][oi]);
						v[hci].append(vos[nci][oi]);
					}
				}
			}
		}
	}
}
