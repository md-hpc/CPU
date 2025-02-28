void simulation::velocity_update(int hci, int tid) {
	voxel hcv = voxelof(hci);

	int nr = s.r[hci].x.size();
	


	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				if (di < 0 || di == 0 && dj < 0 || di == 0 && dj == 0 && dk < 0)
					continue;

				int nci = cell(hcv.i + di, hcv.j + dj, hcv.k + dk);				
				int nn = s.r[nci].size();
				
				for (int ri = 0; ri < nr; ri++) {
					fvec vxa, vya, vza;
					vxa = vya = vza = FVZ;

					rx = s.r[hci].x[ri];
					ry = s.r[hci].y[ri];
					rz = s.r[hci].z[ri];
					for (int p = 0; p < VSIZE; p++) {
						for (int ni = 0; ni < nn; ni++) {
							// nci >= hci from conditional before
							//
							// ni > ri (&& nci == hci) should be omitted to leverag n3l
							// ni == ri (&& nci == hci) requires special handling (no n3l)
							if (nci == hci && ni >= ri) 
								continue;
							fvec vx, vy, vz, vf, r;
							vx = vy = vz = FVZ;

							nx = s.r[nci].x[ni];
							ny = s.r[nci].y[ni];
							nz = s.r[nci].z[ni];
							
							vx = submv(rx, nx);
							vy = submv(ry, ny);
							vz = submv(rz, nz);
							
							r = sqrtv(vx * vx + vy * vy + vz * vz);
							vf = ljv(r);

							vx = vx * vf * vdt / r;
							vy = vy * vf * vdt / r;
							vz = vz * vf * vdt / r;

							vxa += vx;
							vya += vy;
							vza += vz;

							vx *= -1;
							vy *= -1;
							vz *= -1;

							s.v[nci].x[ni] += vx;
							s.v[nci].y[ni] += vy;
							s.v[nci].z[ni] += vz;
						}
				
						permute(&rx, 1);
						permute(&ry, 1);
						permute(&rz, 1);

						permute(&vxa, 1);
						permute(&vya, 1);
						permute(&vza, 1);
					}
					s.v[hci].x[ri] += vxa;
					s.v[hci].y[ri] += vya;
					s.v[hci].y[ri] += vza;
				}
			}
		}
	}

	// handle ri == ni && hci == nci for this cell. Do not apply n3l
	for (int i = 0; i < nr; i++) {
		fvec rx, ry, rz;
		fvec nx, ny, nz;
		fvec vf;
		fvec rva, rya, rza;
		fvec vx, vy, vz;

		rx = s.r[hci].x[i];
		ry = s.r[hci].y[i];
		rz = s.r[hci].z[i];

		nx = s.r[hci].x[i];
		ny = s.r[hci].y[i];
		nz = s.r[hci].z[i];

		nx = permute(nx, 1);
		ny = permute(ny, 1);
		nz = permute(nz, 1);

		for (int p = 0; p < 7; p++) {
			permute(&nx,1);
			permute(&ny,1);
			permute(&nz,1);
			
			vx = submv(rx, nx);
			vy = submv(ry, ny);
			vz = submv(rz, nz);
			
			vf = vx * vx + vy * vy + vz * vz;
			vf = ljv(vf);

			vx *= vf;
			vy *= vf;
			vz *= vz;
			
			vxa += vx;
			vya += vy;
			vza += vz;
		}

		s.v[hci].x[i] += vxa;
		s.v[hci].y[i] += vya;
		s.v[hci].z[i] += vza;
	}
}

void simulation::position_update(int hci, int tid) {
	const float dt = DT;

	int np = s.r[hci].x.size();
	int cur = 0;
	int ocur = 0;

	s.vos[hci].resize(0);
	s.ros[hci].resize(0);
	s.cios[hci].resize(0);

	for (int pi = 0; pi < np; p++) {
		fvec_vector_vec *r, *v;
		r = &s.r[hci];
		v = &s.v[hci]
		r->x[pi] += v->x[pi] * dt;
		r->y[pi] += v->y[pi] * dt;
		r->z[pi] += v->z[pi] * dt;
		
		r->x[pi] = apbcfv(r->x[pi]);
		r->y[pi] = apbcfv(r->y[pi]);
		r->z[pi] = apbcfv(r->z[pi]);

		ipack cells = { .v=cellv(r->x, r->y, r->z) };

		// I think this could be optimized if we tried to do as many aligned vector loads as possible when
		// all of cells == hci
		for (int i = pi * VSIZE; i < (pi + 1) * VSIZE; i++) {
			if (cells.d[i%VSIZE] == hci) {
				r->set(r->get(i),cur);
				v->set(v->get(i),cur);
				cur++;
			} else {
				s.ros[hci].push_back(r->get(i));
				s.vos[hci].push_back(v->get(i));
				s.cios[hci].push_back(cells.d[i%VSIZE]);
			}
		}
	}
	v->resize(cur);
	r->resize(cur);
}

void simulation::cell_update(int hci, int tid) {
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
