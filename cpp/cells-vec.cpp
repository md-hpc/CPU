#include <stdio.h>
#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

using namespace std;

ALGO = ALGO_CELLS_VEC;

typedef int ivec __attribute__ ((vector_size(8)));
typedef float fvec __attribute__((vector_size(8)));
#define VBYTES sizeof(fvec)
#define VSIZE (sizeof(fvec)/sizeof(float))

typedef union {
	fvec v;
	float d[VSIZE];
} pack;

#define VI(v,i) (((pack*)&(v))->d[i])
#define FVZ {0,0,0,0,0,0,0,0}

class fvec_vector_vec {
public:
	vector<fvec> x;
	vector<fvec> y;
	vector<fvec> z;

	fvec_vector_vec() : vi(0) {}

	void append(vec r) {
		fvec v;	
		if (vi == 0) {
			x.push_back(v);
			y.push_back(v);
			z.push_back(v);
		}

		int end = x.size() - 1;
		VI(x[end],vi) = r.x;
		VI(y[end],vi) = r.y;
		VI(z[end],vi) = r.z;

		vi = (vi+1)%VSIZE;
	}

	void dump(vector<vec> v) {
		int n = x.size();
		int i;
		
		int res = vi != 0 ? VSIZE - vi : 0;

		int ai, si;
		for (i = 0; i < VSIZE * n + res; i++) {
			ai = i / VSIZE;
			si = i % VSIZE;
			v.push_back(
				vec(
					VI(x[ai],si), 
					VI(y[ai],si),
					VI(z[ai],si)
				)
			);
		}
	}

private:
	int vi;
}

class simulation {
public:
	vector<fvec_vector_vec> r;
	vector<fvec_vector_vec> v;
}

class voxel_index {
public:
	voxel_index(int i, int j, int k) : i(i), j(j), k(k) {}

	int x;
	int y;
	int z;
}

simulation s;
vector<voxel_index> offsets;


int sign(int i) {
	return i < 0 ? -1 : i == 0 ? 0 : 1;
}

float i2d(int i) {
	return R * (((float)i) - sign(i) * .5)
}

int main(int argc, char **argv) {
	simulation s;
	int np;
	particle *p;

	parse_cli(argc, argv);

	init_particles(particles);

	for (int di = -1; di <= 1; di++) {
		for (int dj = -1; dj <= 1; dj++) {
			for (int dk = -1; dk <= 1; dk++) {
				if (di==0 && dj==0 && dk==0)
					continue;
				offsets.push_back((voxel_index(di,dj,dk));
			}
		}
	}

	np = particles.size();
	for (int pi = 0; pi < np; pi++) {
		p = &particles[pi];
		s.r[p->cell].append(p->r);
		s.v[p->cell].append(p->v);
	}

	for (t = 0; t < TIMESTEPS) {
		thread(velocity_update, N_CELL);
		thread(motion_update, N_CELL);
		thread(cell_update, N_CELL);
	}
}

void permute(fvec *x, int n) {
	pack tmp, *dst;
	tmp.v = *x;
	dst = (pack*) x;
	for (int i = 0; i < VISZE; i++) {
		dst->d[(i+n)%VSIZE] = tmp.d[i];
	}
}

void velocity_update(int hci, int tid) {
	int cci[3], i, j, k;

	fvec vdt = {DT, DT, DT, DT, DT, DT, DT, DT};

	linear_idx(cci, hi);
	i = cci[0];
	j = cci[1];
	k = cci[2];

	int nr = s.r[hci].x.size();
	int no = offsets.size();
	
	for (int oi = 0; oi < no; oi++) {
		int di = offset.i;
		int dj = offset.j;
		int dk = offset.k;
		nci = linear_idx(i+di, j+dj, k+dk);
		
		if (nci < hci)
			continue;

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
					// nci > hci should be omitted to leverage n3l
					// ni > ri (&& nci == hci) should be omitted to leverag n3l
					// ni == ri (&& nci == hci) requires special handling (no n3l)
					if (nci > hci || ni >= ri) 
						continue;
					fvec vx, vy, vz, vf, r;
					vx = vy = vz = FVZ;

					nx = s.r[nci].x[ni];
					ny = s.r[nci].y[ni];
					nz = s.r[nci].z[ni];
					
					// need to do something to verify that this actually works
					vx = submv(rx, nx);
					vy = submv(ry, ny);
					vz = submv(rz, nz);
					
					r = sqrt(vx * vx + vy * vy + vz * vz);
					vf = LJ(r);

					vx = vx * vf * vdt / r;
					vy = vy * vf * vdt / r;
					vz = vz * vf * vdt / r;

					vxa += vx;
					vya += vy;
					vza += vz;

					vx *= -1;
					vy *= -1;
					vz *= -1;

					VLD(s.v[nci].x[ni]) += vx;
					VLD(s.v[nci].y[ni]) += vy;
					VLD(s.v[nci].z[ni]) += vz;
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

			vx = subm(rx, nx);
			vy = subm(ry, ny);
			vz = subm(rz, nz);
			r = sqrt(vx * vx + vy * vy + vz * vz);
			vf = LJ(r) * vdt / r;

			vx *= vf;
			vy *= vf;
			vz *= vz;
			
			vxa += vx;
			vya += vy;
			vza += vz;
		}
	}
}

void motion_update(int hci, int tid) {
	int np = s.r[hci].x.size();

	for (int pi = 0; pi < np; p++) {

	}
}
