// 3D World - Erosion simulation algorithm
// by Frank Gennari
// 2/10/18

#include "3DWorld.h"
#include "mesh.h"


extern float erode_amount, water_plane_z;


// see http://ranmantaru.com/blog/2011/10/08/water-erosion-on-heightmap-terrain/
void apply_erosion(float *heightmap, int xsize, int ysize, float min_zval, unsigned num_iters) {

	if (num_iters == 0 || erode_amount <= 0.0) return; // erosion disabled
	RESET_TIME;
	// Kq and minSlope are for soil carry capacity.
	// Kw is water evaporation speed.
	// Kr is erosion speed (how fast the soil is removed).
	// Kd is deposition speed (how fast the extra sediment is dropped).
	// Ki is direction inertia. Higher values make channel turns smoother.
	// g is gravity that accelerates the flows.
	float const Kq=10, Kw=0.001f, Kr=0.9f, Kd=0.02f, Ki=0.1f, minSlope=0.05f, g=20, Kg=g*2;
	int const PAD(4), NX(xsize+2*PAD), NY(ysize+2*PAD);
	unsigned const MAX_PATH_LEN(4*NX*NY);
	vector<vector2d> erosion(NX*NY, vector2d(0.0, 0.0));
	vector<float> mh_padded(NX*NY);

	// pad mesh by 1 unit on each side to create a buffer of trash around the edges that can be discarded
	for (int y = 0; y < NY; ++y) {
		int const offset(max(min(y-PAD, ysize-1), 0)*xsize);

		for (int x = 0; x < NX; ++x) {
			mh_padded[y*NX + x] = heightmap[max(min(x-PAD, xsize-1), 0) + offset];
		}
	}

#define HMAP_INDEX(x, y) (NX*max(min(y, NY-1), 0) + max(min(x, NX-1), 0))
#define HMAP(x, y) mh_padded[HMAP_INDEX(x, y)]

#define DEPOSIT_AT(X, Z, W) { \
	float const delta = ds*erode_amount*(W); \
	unsigned const ix(HMAP_INDEX((X), (Z))); \
	erosion[ix].y += delta; \
	if (!(X < 0 || Z < 0 || X >= NX || Z >= NY)) {mh_padded[ix] += delta;} \
}

#define DEPOSIT(H) \
	DEPOSIT_AT(xi  , zi  , (1-xf)*(1-zf)) \
	DEPOSIT_AT(xi+1, zi  ,    xf *(1-zf)) \
	DEPOSIT_AT(xi  , zi+1, (1-xf)*   zf ) \
	DEPOSIT_AT(xi+1, zi+1,    xf *   zf ) \
	(H)+=ds;

#define ERODE(X, Z, W) { \
	float const delta=ds*erode_amount*(W); \
	unsigned const ix(HMAP_INDEX((X), (Z))); \
	mh_padded[ix]-=delta; \
	vector2d &e=erosion[ix]; \
	float r=e.x, d=e.y; \
	if (delta<=d) {d-=delta;} else {r+=delta-d; d=0;} \
	e.x=r; e.y=d; \
}

#pragma omp parallel for schedule(dynamic,1)
	for (int iter=0; iter < (int)num_iters; ++iter) {
		rand_gen_t rgen;
		rgen.set_state(iter+11, 79*iter+121);
		int xi = PAD + (rgen.rand()%xsize);
		int zi = PAD + (rgen.rand()%ysize);
		float xp=xi, zp=zi, xf=0, zf=0, s=0, v=0, w=1, dx=0, dz=0;
		float h=HMAP(xi, zi), h00=h, h10=HMAP(xi+1, zi), h01=HMAP(xi, zi+1), h11=HMAP(xi+1, zi+1);

		unsigned numMoves=0;
		for (; numMoves<MAX_PATH_LEN; ++numMoves) {
			// calc gradient
			float gx=h00+h01-h10-h11, gz=h00+h10-h01-h11;
			// calc next pos
			dx=(dx-gx)*Ki+gx;
			dz=(dz-gz)*Ki+gz;

			float dl=sqrtf(dx*dx+dz*dz);
			if (dl<=FLT_EPSILON) { // pick random dir
				float a=rgen.rand_float()*TWO_PI;
				dx=cosf(a); dz=sinf(a);
			}
			else {
				dx/=dl; dz/=dl;
			}
			float nxp=xp+dx, nzp=zp+dz;
			// sample next height
			int nxi=floor(nxp), nzi=floor(nzp);
			float nxf=nxp-nxi, nzf=nzp-nzi;
			float nh00=HMAP(nxi, nzi), nh10=HMAP(nxi+1, nzi), nh01=HMAP(nxi, nzi+1), nh11=HMAP(nxi+1, nzi+1);
			float nh=(nh00*(1-nxf)+nh10*nxf)*(1-nzf)+(nh01*(1-nxf)+nh11*nxf)*nzf;
			// adjust by HALF_DXY = average mesh texel size - this is river depth
			if (max(max(nh00, nh10), max(nh01, nh11)) < water_plane_z - HALF_DXY) break; // reached ocean water, stop and ignore sediment

			// if higher than current, try to deposit sediment up to neighbour height
			bool const outside(xi < 0 || zi < 0 || xi >= NX || zi >= NY);
			if (nh>=h || outside) {
				float ds=(nh-h)+0.001f;

				if (ds>=s || outside) {
					ds=s;
					DEPOSIT(h) // deposit all sediment
					s=0;
					break; // stop
				}
				DEPOSIT(h)
				s-=ds;
				v=0;
			}
			// compute transport capacity
			float dh=h-nh;
			float slope=dh;
			//float slope=dh/sqrtf(dh*dh+1);
			float q=max(slope, minSlope)*v*w*Kq;

			// deposit/erode (don't erode more than dh)
			float ds=s-q;
			if (ds>=0) { // deposit
				ds*=Kd;
				//ds=minval(ds, 1.0f);
				DEPOSIT(dh)
				s-=ds;
			}
			else { // erode
				ds*=-Kr;
				ds=min(ds, dh*0.99f);
				ds*=((get_bare_ls_tid(nh) == ROCK_TEX) ? 0.5 : 2.0); // rock erodes slower than dirt/sand

				for (int z=zi-1; z<=zi+2; ++z) {
					float zo=z-zp, zo2=zo*zo;

					for (int x=xi-1; x<=xi+2; ++x) {
						float xo=x-xp;
						float w=1-(xo*xo+zo2)*0.25f;
						if (w<=0) continue;
						w*=0.1591549430918953f;
						ERODE(x, z, w)
					}
				}
				dh-=ds;
				s+=ds;
			}
			// move to the neighbor
			v=sqrtf(v*v+Kg*dh);
			w*=1-Kw;
			xp=nxp; zp=nzp; xi=nxi; zi=nzi; xf=nxf; zf=nzf;
			h=nh; h00=nh00; h10=nh10; h01=nh01; h11=nh11;
		} // for numMoves
		if (numMoves>=MAX_PATH_LEN) {cout << "droplet path is too long: " << iter << endl;}
	} // for iter

	// remove padding and clamp to min_zval
	for (int y = 0; y < ysize; ++y) {
		for (int x = 0; x < xsize; ++x) {
			heightmap[y*xsize + x] = max(min_zval, mh_padded[(y+PAD)*NX + x+PAD]);
		}
	}
	PRINT_TIME("Erosion");
}

