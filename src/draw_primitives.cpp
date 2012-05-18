// 3D World - Shape primitive drawing
// by Frank Gennari
// 3/6/06
#include "3DWorld.h"
#include "subdiv.h"
#include "upsurface.h"
#include "gl_ext_arb.h"


// predefined display lists
enum {DLIST_CYLIN=0, DLIST_CYLIN_T, DLIST_CONE, DLIST_CONE_T, NUM_RES_DLIST};

bool const USE_SPHERE_DLIST = 1;
bool const USE_CYLIN_DLIST  = 1;

unsigned const NUM_SPH_DLIST(4*N_SPHERE_DIV), NUM_CYLIN_DLIST(N_CYL_SIDES);

unsigned predef_dlists[NUM_RES_DLIST]   = {0};
unsigned sphere_dlists[NUM_SPH_DLIST]   = {0};
unsigned cylin_dlists [NUM_CYLIN_DLIST] = {0};
vector_point_norm cylinder_vpn;


extern int display_mode, draw_model;
extern GLUquadricObj* quadric;


#define NDIV_SCALE(val) unsigned(val*ndiv + 0.5)


// ******************** CURVE GENERATORS ********************


void get_ortho_vectors(vector3d const &v12, vector3d *vab, int force_dim) {

	assert(vab != NULL);
	vector3d vtest(v12);
	unsigned dim(0);

	if (force_dim >= 0) {
		assert(force_dim < 3);
		dim = force_dim;
	}
	else {
		dim = get_min_dim(v12);
	}
	vtest[dim] += 0.5;
	cross_product(vtest, v12,    vab[0]); // vab[0] is orthogonal to v12
	cross_product(v12,   vab[0], vab[1]); // vab[1] is orthogonal to v12 and vab[0]
	for (unsigned i = 0; i < 2; ++i) vab[i].normalize();
}


// create class sd_cylin_d?
// perturb_map usage is untested
vector_point_norm const &gen_cylinder_data(point const ce[2], float radius1, float radius2, unsigned ndiv, vector3d &v12,
										   float const *const perturb_map, float s_beg, float s_end, int force_dim)
{
	assert(ndiv > 0 && ndiv < 100000);
	v12 = (ce[1] - ce[0]).get_norm();
	float const r[2] = {radius1, radius2};
	float const cs_scale(TWO_PI/(float)ndiv);
	vector3d vab[2];
	get_ortho_vectors(v12, vab, force_dim);
	vector_point_norm &vpn(cylinder_vpn);
	unsigned s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end));
	if (s1 == ndiv) s0 = 0; // make wraparound correct
	s1 = min(ndiv, s1+1);   // allow for sn
	vpn.p.resize(2*ndiv);
	vpn.n.resize(ndiv);
	float const css(TWO_PI/(float)ndiv), sin_ds(sin(css)), cos_ds(cos(css));
	float sin_s((s0 == 0) ? 0.0 : sin(s0*css)), cos_s((s0 == 0) ? 1.0 : cos(s0*css));
	// sin(x + y) = sin(x)*cos(y) + cos(x)*sin(y)
	// cos(x + y) = cos(x)*cos(y) - sin(x)*sin(y)

	for (unsigned S = s0; S < s1; ++S) { // build points table
		float const s(sin_s), c(cos_s);

		for (unsigned i = 0; i < 2; ++i) {
			vpn.p[(S<<1)+i] = ce[i] + (vab[0]*(r[i]*s)) + (vab[1]*(r[i]*c));
		}
		sin_s = s*cos_ds + c*sin_ds;
		cos_s = c*cos_ds - s*sin_ds;
	}
	float nmag_inv(0.0);
	bool const npt_r2(radius1 < radius2); // determine normal from longest edge for highest accuracy (required for when r1 == 0.0)

	for (unsigned S = s0; S < s1; ++S) { // build normals table
		vector3d const v1(vpn.p[(((S+1)%ndiv)<<1)+npt_r2], vpn.p[(S<<1)+npt_r2]), v2(vpn.p[(S<<1)+1], vpn.p[S<<1]);
		cross_product(v2, v1, vpn.n[S]);
		if (S == s0) nmag_inv = 1.0/vpn.n[S].mag(); // first one (should not have a divide by zero)
		vpn.n[S] *= nmag_inv; //norms[S].normalize();
	}
	if (perturb_map != NULL) { // add in perturbations
		float const ravg(0.5*(r[0] + r[1]));
		float const pscale[2] = {r[0]/ravg, r[1]/ravg};

		for (unsigned S = s0; S < s1; ++S) {
			for (unsigned i = 0; i < 2; ++i) {
				vpn.p[(S<<1)+i] += vpn.n[S]*(pscale[i]*perturb_map[S]);
			}
		}
		// update normals?
	}
	return vpn;
}


void sd_sphere_d::gen_points_norms(float s_beg, float s_end, float t_beg, float t_end) {

	unsigned const ndiv(spn.ndiv);
	assert(ndiv > 0 && ndiv < 100000);
	float const cs_scale(PI/(float)ndiv), cs_scale2(2.0*cs_scale), sin_dt(sin(cs_scale)), cos_dt(cos(cs_scale));
	static sphere_point_norm temp_spn;

	if (temp_spn.points == NULL || ndiv > temp_spn.ndiv) { // allocate all memory
		if (temp_spn.points != NULL) temp_spn.free(); // already allocated at least once
		temp_spn.alloc(ndiv);
		unsigned const stride(ndiv+1);
	}
	else {
		temp_spn.set_pointer_stride(ndiv);
	}
	unsigned s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end)), t0(NDIV_SCALE(t_beg)), t1(NDIV_SCALE(t_end));
	if (s1 == ndiv) s0 = 0; // make wraparound correct
	s1 = min(ndiv, s1+1);   // allow for sn
	t1 = min(ndiv, t1+1);   // allow for tn
	float sin_t0((t0 == 0) ? 0.0 : sin(t0*cs_scale)), cos_t0((t0 == 0) ? 1.0 : cos(t0*cs_scale));

	for (unsigned s = s0; s < s1; ++s) { // build points and normals table
		float const theta(s*cs_scale2), tvc(cos(theta)), tvs(sin(theta)); // theta1, theta2
		float sin_t(sin_t0), cos_t(cos_t0);
		unsigned const soff((ndiv+1)*s);

		for (unsigned t = t0; t <= t1; ++t) { // Note: x and y are swapped because theta is out of phase by 90 degrees to match tex coords
			point &pt(temp_spn.points[s][t]);
			float const sv(sin_t), cv(cos_t);
			pt.assign(sv*tvs, sv*tvc, cv); // R*sin(phi)*sin(theta), R*sin(phi)*cos(theta), R*cos(phi)
			sin_t = sv*cos_dt + cv*sin_dt;
			cos_t = cv*cos_dt - sv*sin_dt;
			temp_spn.norms[s][t] = pt; // Note: perturb_map does not affect normals until later
			pt   *= radius;
			pt   += pos;
			if (perturb_map) pt += temp_spn.norms[s][t]*perturb_map[t+soff];
			if (surf)        pt += temp_spn.norms[s][t]*surf->get_height_at(pt);
		}
	}
	if (perturb_map || surf) { // recalculate vertex/surface normals
		for (unsigned s = s0; s < s1; ++s) {
			for (unsigned t = max(1U, t0); t <= min(t1, ndiv-1); ++t) { // skip the poles
				point const p(temp_spn.points[s][t]);
				vector3d n(zero_vector);

				for (unsigned d = 0; d < 4; ++d) { // s+,t+  t-,s+  s-,t-  t+,s-
					unsigned const si((d==0||d==1) ? (s+1)%ndiv : (s+ndiv-1)%ndiv);
					unsigned const ti((d==0||d==3) ? (t+1)      : (t-1));
					vector3d const nst(cross_product((temp_spn.points[si][t] - p), (temp_spn.points[s][ti] - p)));
					float const nmag(nst.mag());
					if (nmag > TOLERANCE) n += nst*(((d&1) ? -0.25 : 0.25)/nmag);
				}
				temp_spn.norms[s][t] = n;
			}
		}
	}
	spn.points = temp_spn.points;
	spn.norms  = temp_spn.norms;
}


float sd_sphere_d::get_rmax() const { // could calculate this during gen_points_norms

	assert(spn.points);
	float rmax_sq(0.0);

	for (unsigned y = 0; y < spn.ndiv; ++y) {
		for (unsigned x = 0; x <= spn.ndiv; ++x) {
			rmax_sq = max(rmax_sq, p2p_dist_sq(pos, spn.points[y][x]));
		}
	}
	return sqrt(rmax_sq);
}


// *** sd_sphere_d ***


void sd_sphere_d::set_data(point const &p, float r, int n, float const *pm, float dp, upsurface const *const s) {

	pos         = p;
	radius      = r;
	def_pert    = dp;
	perturb_map = pm;
	surf        = s;
	spn.ndiv    = n;
	assert(radius > 0.0);
}


void sd_sphere_d::make_local_copy() {

	assert(!local);
	local = 1;
	point *pdata(spn.points[0]); // make a temporary copy
	spn.alloc(spn.ndiv);
	assert(sizeof(point) == sizeof(vector3d));
	memcpy(spn.points[0], pdata, 2*(spn.ndiv*(spn.ndiv+1))*sizeof(point));
}


void sd_sphere_d::free_local_data() {

	assert(local);
	local = 0;
	spn.free();
}


// *** sphere_point_norm ***


void sphere_point_norm::alloc(unsigned ndiv_) {

	ndiv = ndiv_;
	matrix_base_alloc_2d(points, (ndiv+1), 2*ndiv);
	set_pointer_stride(ndiv);
}


void sphere_point_norm::set_pointer_stride(unsigned ndiv_) {

	ndiv     = ndiv_;
	norms    = points    + ndiv; // assumes point and vector3d are the same
	norms[0] = points[0] + ndiv*(ndiv+1); // num_alloc
	matrix_ptr_fill_2d(points, (ndiv+1), ndiv);
	matrix_ptr_fill_2d(norms,  (ndiv+1), ndiv);
}


void sphere_point_norm::free() {

	assert(points && norms && points[0] && norms[0]);
	matrix_delete_2d(points);
	norms = NULL;
}


// ******************** CYLINDER ********************


void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals) {

	if (invert_normals) gluQuadricOrientation(quadric, GLU_INSIDE);
	gluDisk(quadric, r_inner, r_outer, ndiv, 1);
	if (invert_normals) gluQuadricOrientation(quadric, GLU_OUTSIDE);
}


void draw_circle_normal_at_z(float z, float r_inner, float r_outer, int ndiv, int invert_normals) {

	glPushMatrix();
	glTranslatef(0.0, 0.0, z);
	draw_circle_normal(r_inner, r_outer, ndiv, invert_normals);
	glPopMatrix();
}


void draw_cylinder(float length, float radius1, float radius2, int ndiv, int nstacks,
				   bool draw_ends, bool first_end_only, bool last_end_only)
{ // length can be negative
	assert(quadric && ndiv > 0 && nstacks > 0);
	//set_fill_mode(); // too slow? not needed?
	//if (nstacks == 1) draw_cylin_fast(radius1, radius2, length, ndiv, 0, 1); // tex coords?
	gluCylinder(quadric, radius1, radius2, length, ndiv, nstacks); // draw quad/triangle if small?
	
	if (draw_ends) {
		if (!last_end_only && radius1 > 0.0) draw_circle_normal(0.0, radius1, ndiv, 1);
		
		if (!first_end_only && radius2 > 0.0) {
			glPushMatrix();
			glTranslatef(0.0, 0.0, length);
			draw_circle_normal(0.0, radius2, ndiv, 0);
			glPopMatrix();
		}
	}
}


void draw_cylinder_nstacks(float len, float r1, float r2, int ndiv, int nstacks, bool texture) { // no draw_ends

	if (nstacks == 1) {
		draw_cylin_fast(r1, r2, len, ndiv, texture, 1);
	}
	else {
		draw_cylinder(len, r1, r2, ndiv, nstacks, 0);
	}
}


void draw_cylinder(point const &p1, float length, float radius1, float radius2, int ndiv, int nstacks, bool draw_ends) {

	glPushMatrix();
	translate_to(p1);
	draw_cylinder(length, radius1, radius2, ndiv, nstacks, draw_ends);
	glPopMatrix();
}


void xform_cylinder(point const &p1, point const &p2, vector3d const &scale) {

	translate_to(p1);
	vector3d const v2z(p2, p1);
	rotate_into_plus_z(v2z);
	
	if (scale != zero_vector) {
		vector3d s(scale);
		rotate_vector3d_by_vr(v2z, plus_z, s);
		rotate_from_v2v(s, vector3d(1.0, 0.0, 0.0)); // rotate into +x
		glScalef(s.mag(), 1.0, 1.0);
	}
}


void draw_rotated_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, int nstacks,
						   bool draw_ends, vector3d const &scale)
{
	glPushMatrix();
	xform_cylinder(p1, p2, scale);
	draw_cylinder(all_zeros, p2p_dist(p1, p2), radius1, radius2, ndiv, nstacks, draw_ends);
	glPopMatrix();
}


// draw_sides_ends: 0 = draw sides only, 1 = draw sides and ends, 2 = draw ends only
void draw_fast_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, bool texture,
						int draw_sides_ends, float const *const perturb_map, bool const *const render_map,
						float const *const exp_map, point const *const pt_shift, float expand, float s_beg, float s_end)
{
	assert(radius1 > 0.0 || radius2 > 0.0);
	bool const use_quads(render_map || pt_shift || exp_map || expand != 0.0);
	point const ce[2] = {p1, p2};
	float const ndiv_inv(1.0/ndiv);
	vector3d v12; // (ce[1] - ce[0]).get_norm()
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, ndiv, v12, perturb_map, s_beg, s_end));
	if (expand != 0.0) expand *= 0.5; // 1/2, for normalization
	unsigned const s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end));

	if (draw_sides_ends == 2) {
		// draw ends only - nothing to do here
	}
	else if (use_quads) {
		glBegin(GL_QUADS);
		
		for (unsigned S = s0; S < s1; ++S) {
			if (render_map && !render_map[S]) continue;
			float const exp(expand + (exp_map ? exp_map[S] : 0.0));
			unsigned const s[2] = {S, (S+ndiv-1)%ndiv};

			for (unsigned i = 0; i < 4; ++i) {
				unsigned const ss(s[i>>1]);
				if (texture) glTexCoord2f((1.0 - ss*ndiv_inv), float(i==1||i==2));
				point pt(vpn.p[(ss<<1) + (i==1||i==2)]);
				if (exp != 0.0) pt += (vpn.n[s[0]] + vpn.n[s[1]])*exp; // average the normals
				if (pt_shift) pt += pt_shift[S];
				vpn.n[ss].do_glNormal();
				pt.do_glVertex();
			}
		}
		glEnd();
	}
	else if (radius2 == 0.0) { // cone
		glBegin(GL_TRIANGLES);

		for (unsigned s = s0; s < s1; ++s) {
			unsigned const sp((s+ndiv-1)%ndiv), sn((s+1)%ndiv);
			//if (texture) glTexCoord2f((1.0 - (s+0.5)*ndiv_inv), 1.0); // small discontinuities at every position
			if (texture) glTexCoord2f(0.5, 1.0); // one big discontinuity at one position
			vpn.n[s].do_glNormal();
			vpn.p[(s <<1)+1].do_glVertex();
			if (texture) glTexCoord2f((1.0 - (s+0.0)*ndiv_inv), 0.0);
			(vpn.n[s] + vpn.n[sp]).do_glNormal(); // normalize?
			vpn.p[(s <<1)+0].do_glVertex();
			if (texture) glTexCoord2f((1.0 - (s+1.0)*ndiv_inv), 0.0);
			(vpn.n[s] + vpn.n[sn]).do_glNormal(); // normalize?
			vpn.p[(sn<<1)+0].do_glVertex();
		}
		glEnd();
	}
	else {
		glBegin(GL_QUAD_STRIP);

		for (unsigned S = s0; S <= s1; ++S) {
			unsigned const s(S%ndiv);
			(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]).do_glNormal(); // normalize?

			for (unsigned i = 0; i < 2; ++i) {
				if (texture) glTexCoord2f((1.0 - S*ndiv_inv), float(i));
				vpn.p[(s<<1)+i].do_glVertex();
			}
		}
		glEnd();
	}
	if (draw_sides_ends != 0) {
		assert(!render_map && !exp_map && !pt_shift && expand == 0.0 && s0 == 0 && s1 == ndiv); // not yet supported
		float const r[2] = {radius1, radius2};

		for (unsigned i = 0; i < 2; ++i) {
			if (r[i] == 0.0) continue;
			(v12*(i ? 1.0 : -1.0)).do_glNormal();
			glBegin(GL_TRIANGLE_FAN);
			if (texture) glTexCoord2f(0.5, 0.5);
			ce[i].do_glVertex();

			for (unsigned S = s0; S <= s1; ++S) { // ndiv can change
				unsigned const ss(S%ndiv), s(i ? (ndiv - ss - 1) : ss);
				
				if (texture) { // inefficient
					float const theta(TWO_PI*s/ndiv);
					glTexCoord2f(0.5*(1.0 + sinf(theta)), (0.5*(1.0 + cosf(theta))));
				}
				vpn.p[(s<<1)+i].do_glVertex();
			}
			glEnd();
		}
	}
}


void draw_cylindrical_section(point const &pos, float length, float r_inner, float r_outer, int ndiv, bool texture) {

	assert(r_outer > 0.0 && r_inner >= 0.0 && length >= 0.0 && ndiv > 0 && r_outer >= r_inner);
	glPushMatrix();
	translate_to(pos);
	draw_cylin_fast(r_outer, r_outer, length, ndiv, texture, 1);

	if (r_inner != r_outer) {
		if (r_inner != 0.0) draw_cylin_fast(r_inner, r_inner, length, ndiv, texture, 1);
		draw_circle_normal(r_inner, r_outer, ndiv, 1); // no texture, could make loops > 1
		glTranslatef(0.0, 0.0, length);
		draw_circle_normal(r_inner, r_outer, ndiv, 0); // no texture, could make loops > 1
	}
	glPopMatrix();
}


void draw_trunc_cone(point pos, vector3d v1, float length, float radius, float radius2, bool is_camera) {

	int const nsides(N_CYL_SIDES >> (is_camera ? 0 : 1));
	glPushMatrix();
	translate_to(pos);
	rotate_by_vector(v1, 180.0);
	gluCylinder(quadric, radius, radius2, length, nsides, 1);
	if (is_camera || sphere_in_camera_view(pos, radius, 0)) gluDisk(quadric, 0.0, radius, nsides, 1);
	glPopMatrix();
	pos += v1*length; // destroys pos

	if (is_camera || sphere_in_camera_view(pos, radius2, 0)) {
		glPushMatrix();
		translate_to(pos);
		rotate_by_vector(v1, 180.0);
		gluDisk(quadric, 0.0, radius2, nsides, 1);
		glPopMatrix();
	}
}


// ******************** SPHERE ********************


void draw_sphere_at(point const &pos, float radius, int ndiv) {

	assert(quadric && radius > 0.0);
	set_fill_mode(); // might remove later
	bool const nonzero(pos != all_zeros);
	if (nonzero) glPushMatrix();
	if (nonzero) translate_to(pos);
	gluSphere(quadric, radius, ndiv, ndiv);
	if (nonzero) glPopMatrix();
}


void draw_sphere_at_tc(point const &pos, float radius, int ndiv, bool texture, bool cull) {

	if (texture) gluQuadricTexture(quadric, GL_TRUE);
	if (cull)    glEnable(GL_CULL_FACE);
	draw_sphere_at(pos, radius, ndiv);
	if (cull)    glDisable(GL_CULL_FACE);
	if (texture) gluQuadricTexture(quadric, GL_FALSE);
}


// back face culling requires the sphere to be untransformed (or the vertices to be per-transformed)
void sd_sphere_d::draw_subdiv_sphere(point const &vfrom, int texture, bool disable_bfc, bool const *const render_map,
									 float const *const exp_map, point const *const pt_shift, float expand,
									 float s_beg, float s_end, float t_beg, float t_end, unsigned sv1) const
{
	//draw_sphere_at_tc(pos, radius, spn.ndiv, 1, 0);
	assert(!render_map || disable_bfc);
	unsigned const ndiv(spn.ndiv);
	assert(ndiv > 0 && sv1 > 0);
	float const ndiv_inv(1.0/float(ndiv)), rv(def_pert + radius), rv_sq(rv*rv), tscale(texture);
	float const toler(1.0E-6*radius*radius + rv_sq), dmax(rv + 0.1*radius), dmax_sq(dmax*dmax);
	point **points   = spn.points;
	vector3d **norms = spn.norms;
	bool const use_quads(render_map || pt_shift || exp_map || expand != 0.0);
	if (use_quads) glBegin(GL_QUADS);
	if (expand != 0.0) expand *= 0.25; // 1/4, for normalization
	unsigned const s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end)), t0(NDIV_SCALE(t_beg)), t1(NDIV_SCALE(t_end));
	
	for (unsigned s = s0; s < s1; s += sv1) {
		s = min(s, s1-1);
		unsigned const sn((s+sv1)%ndiv), snt(min((s+sv1), s1));

		if (use_quads) { // use slower quads - no back face culling
			for (unsigned t = t0; t < t1; t += sv1) {
				t = min(t, t1);
				unsigned const ix(s*(ndiv+1)+t);
				if (render_map && !render_map[ix]) continue;
				float const exp(expand + (exp_map ? exp_map[ix] : 0.0));
				unsigned const tn(min(t+sv1, ndiv+1));
				point          pts[4]     = {points[s][t], points[sn][t], points[sn][tn], points[s][tn]};
				vector3d const normals[4] = {norms [s][t], norms [sn][t], norms [sn][tn], norms [s][tn]};
				
				if (exp != 0.0) { // average the normals
					vector3d const quad_norm(normals[0] + normals[1] + normals[2] + normals[3]);
					for (unsigned i = 0; i < 4; ++i) pts[i] += quad_norm*exp;
				}
				for (unsigned i = 0; i < 4; ++i) {
					if (texture) {
						glTexCoord2f(tscale*(1.0f - (((i&1)^(i>>1)) ? snt : s)*ndiv_inv),
                                     tscale*(1.0f - ((i>>1) ? tn : t)*ndiv_inv));
					}
					if (pt_shift) pts[i] += pt_shift[ix];
					normals[i].do_glNormal();
					pts[i].do_glVertex();
				}
			}
		}
		else { // use triangle strips
			glBegin(GL_TRIANGLE_STRIP);

			for (unsigned t = t0; t <= t1; t += sv1) {
				t = min(t, t1);
				point    const pts[2]     = {points[s][t], points[sn][t]};
				vector3d const normals[2] = {norms [s][t], norms [sn][t]};
				bool draw(disable_bfc);

				for (unsigned d = 0; d < 2 && !draw; ++d) {
					float const dp(dot_product_ptv(normals[d], vfrom, pts[d]));
					
					if (dp >= 0.0) {
						draw = 1;
					}
					else if (perturb_map != NULL) { // sort of a hack, not always correct but good enough in most cases
						float const dist_sq(p2p_dist_sq(pos, pts[d]));
						if (dist_sq > toler && (dist_sq > dmax_sq || dp > -0.3*p2p_dist(vfrom, pts[d]))) draw = 1;
					}
				}
				if (draw) {
					for (unsigned i = 0; i < 2; ++i) {
						if (texture) glTexCoord2f(tscale*(1.0f - (i ? snt : s)*ndiv_inv), tscale*(1.0f - t*ndiv_inv));
						normals[i].do_glNormal();
						pts[i].do_glVertex();
					}
				}
			}
			glEnd();
		}
	}
	if (use_quads) glEnd();
}


void sd_sphere_d::draw_ndiv_pow2(unsigned ndiv) const {

	ndiv = max(ndiv, 4U);
	unsigned skip(1);
	for (unsigned n = (spn.ndiv >> 1); ndiv < n; n >>= 1, skip <<= 1) {}
	draw_subdiv_sphere(zero_vector, 1, 1, NULL, NULL, NULL, 0.0, 0.0, 1.0, 0.0, 1.0, skip); // no bfc
}


// non-collision object version
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, point const &vfrom, float const *perturb_map,
						int texture, bool disable_bfc, bool const *const render_map, float const *const exp_map,
						point const *const pt_shift, float expand, float s_beg, float s_end, float t_beg, float t_end)
{
	sd_sphere_d sd(pos, radius, ndiv, perturb_map);
	sd.gen_points_norms(s_beg, s_end, t_beg, t_end);
	sd.draw_subdiv_sphere(vfrom, texture, disable_bfc, render_map, exp_map, pt_shift, expand, s_beg, s_end, t_beg, t_end);
}


void draw_subdiv_sphere(point const &pos, float radius, int ndiv, int texture, bool disable_bfc) {

	draw_subdiv_sphere(pos, radius, ndiv, (disable_bfc ? all_zeros : get_camera_all()), NULL, texture, disable_bfc);
}


void draw_subdiv_sphere_section(point const &pos, float radius, int ndiv, int texture,
								float s_beg, float s_end, float t_beg, float t_end)
{
	draw_subdiv_sphere(pos, radius, ndiv, all_zeros, NULL, texture, 1, NULL, NULL, NULL, 0.0, s_beg, s_end, t_beg, t_end);
}


void rotate_sphere_tex_to_dir(vector3d const &dir) { // dir must be normalized

	glRotatef(atan2(-dir.y, -dir.x)*TO_DEG-90.0, 0.0, 0.0, 1.0); // negate because dir was backwards
	glRotatef(asinf(-dir.z)*TO_DEG, 1.0, 0.0, 0.0);
}


// unused
// less efficient than regular sphere
// textures must tile along all edges for this to look correct
void draw_cube_map_sphere(point const &pos, float radius, int ndiv, bool texture, bool disable_bfc) {

	point pt;
	float const step(1.0/ndiv);
	point const camera(get_camera_all());

	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			pt[n] = (float)j - 0.5;

			for (int s = 0; s < ndiv; ++s) {
				float const va[2] = {(step*(s+1)-0.5), (step*s-0.5)};
				glBegin(GL_QUAD_STRIP);

				for (int t = 0; t <= ndiv; ++t) {
					float s[2];
					pt[d[1]] = s[1] = step*t - 0.5;

					if (!disable_bfc) {
						bool back_facing(1);

						for (unsigned k = 0; k < 2 && back_facing; ++k) {
							point pt2(pt);
							pt2[d[0]] = va[k];
							if (dot_product_ptv(pt2, camera, (pt2.get_norm()*radius + pos)) > 0.0) back_facing = 0;
						}
						if (back_facing) continue;
					}
					for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
						pt[d[0]] = s[0] = va[k^j]; // need to orient the vertices differently for each side
						if (texture) glTexCoord2fv(s);
						vector3d const norm(pt.get_norm());
						norm.do_glNormal();
						(norm*radius + pos).do_glVertex();
					}
				}
				glEnd();
			}
		}
	}
}


// ******************** TORUS ********************


void draw_torus(float ri, float ro, unsigned ndivi, unsigned ndivo, bool do_tex) { // at (0,0,0) in z-plane

	assert(ndivi > 2 && ndivo > 2);
	float const ts(1.0/float(ndivo)), tt(1.0/float(ndivi)), ds(TWO_PI*ts), dt(TWO_PI*tt), cds(cos(ds)), sds(sin(ds));
	static vector<float> sin_cos;
	sin_cos.resize(2*ndivi);

	for (unsigned t = 0; t < ndivi; ++t) {
		float const phi(t*dt);
		sin_cos[(t<<1)+0] = cos(phi);
		sin_cos[(t<<1)+1] = sin(phi);
	}
	for (unsigned s = 0; s < ndivo; ++s) {
		float const theta(s*ds), ct(cos(theta)), st(sin(theta));
		point const pos[2] = {point(ct, st, 0.0), point((ct*cds - st*sds), (st*cds + ct*sds), 0.0)};
		glBegin(GL_QUAD_STRIP);

		for (unsigned t = 0; t <= ndivi; ++t) {
			unsigned const t_((t == ndivi) ? 0 : t);
			float const cp(sin_cos[(t_<<1)+0]), sp(sin_cos[(t_<<1)+1]);

			for (unsigned i = 0; i < 2; ++i) {
				if (do_tex) glTexCoord2f(ts*(s+1-i), tt*t);
				vector3d const delta(point(0.0, 0.0, cp) + pos[1-i]*sp);
				delta.do_glNormal();
				(pos[1-i]*ro + delta*ri).do_glVertex();
			}
		}
		glEnd();
	}
}


// ******************** QUAD/CUBE/POLYGON ********************


void rotate_towards_camera(point const &pos) {

	rotate_into_plus_z((get_camera_pos() - pos));
}


void draw_textured_square(float size, float z, int tid) {

	draw_textured_quad(size, size, z, tid);
}


void draw_textured_square_alpha_test(float size, float z, int tid) {

	GLboolean const blend(glIsEnabled(GL_BLEND));
	if (!blend) enable_blend();
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.001);
	draw_textured_quad(size, size, z, tid);
	glDisable(GL_ALPHA_TEST);
	if (!blend) disable_blend();
}


void draw_flare_no_blend(point const &pos, point const &xlate, float xsize, float ysize) {

	glDepthMask(GL_FALSE);
	point const camera(get_camera_pos());
	select_texture(BLUR_TEX);
	(camera - pos).do_glNormal();
	glBegin(GL_QUADS);
	draw_billboard(xlate, (camera - pos + xlate), up_vector, xsize, ysize);
	glEnd();
	glDepthMask(GL_TRUE);
	glDisable(GL_TEXTURE_2D);
}


void draw_flare(point const &pos, point const &xlate, float xsize, float ysize) {

	enable_blend();
	draw_flare_no_blend(pos, xlate, xsize, ysize);
	disable_blend();
}


void enable_flares(colorRGBA const &color, bool zoomed) { // alpha test?

	glDisable(GL_LIGHTING);
	color.do_glColor();
	glDepthMask(GL_FALSE); // not quite right - prevents flares from interfering with each other but causes later shapes to be drawn on top of the flares
	enable_blend();
	if (draw_model == 0) select_texture(zoomed ? BLUR_CENT_TEX : BLUR_TEX);
	glNormal3f(0.0, 0.0, 1.0);
}


void disable_flares() {

	glDisable(GL_TEXTURE_2D);
	disable_blend();
	glDepthMask(GL_TRUE);
	glEnable(GL_LIGHTING);
}


void draw_textured_quad(float xsize, float ysize, float z, int tid) {

	select_texture(tid);
	plus_z.do_glNormal();
	draw_tquad(xsize, ysize, z, 1);
}


void draw_tquad(float xsize, float ysize, float z, bool texture) {

	glBegin(GL_QUADS);
	draw_one_tquad(-xsize, -ysize, xsize, ysize, z, texture);
	glEnd();
}


void draw_one_tquad(float x1, float y1, float x2, float y2, float z, bool texture, float tx1, float ty1, float tx2, float ty2) {

	if (texture) {glTexCoord2f(tx1, ty1);} glVertex3f(x1, y1, z);
	if (texture) {glTexCoord2f(tx1, ty2);} glVertex3f(x1, y2, z);
	if (texture) {glTexCoord2f(tx2, ty2);} glVertex3f(x2, y2, z);
	if (texture) {glTexCoord2f(tx2, ty1);} glVertex3f(x2, y1, z);
}


void draw_one_mult_tex_quad(float x1, float y1, float x2, float y2, float z, float tx1, float ty1, float tx2, float ty2) {

	multitex_coord2f_a(tx1, ty1); glVertex3f(x1, y1, z);
	multitex_coord2f_a(tx1, ty2); glVertex3f(x1, y2, z);
	multitex_coord2f_a(tx2, ty2); glVertex3f(x2, y2, z);
	multitex_coord2f_a(tx2, ty1); glVertex3f(x2, y1, z);
}


void draw_billboard(point const &pos, point const &viewer, vector3d const &up_dir,
					float xsize, float ysize, float tx1, float ty1, float tx2, float ty2)
{
	vector3d const vdir(viewer - pos);
	vector3d const v1((cross_product(vdir, up_dir).get_norm())*xsize); // what if colinear?
	vector3d const v2((cross_product(vdir, v1    ).get_norm())*ysize);
	glTexCoord2f(tx1, ty1); (pos - v1 - v2).do_glVertex();
	glTexCoord2f(tx1, ty2); (pos - v1 + v2).do_glVertex();
	glTexCoord2f(tx2, ty2); (pos + v1 + v2).do_glVertex();
	glTexCoord2f(tx2, ty1); (pos + v1 - v2).do_glVertex();
}


void draw_line_tquad(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2, bool globalize) {

	int npts(0);
	point const pcenter((p1 + p2)*0.5);
	vector3d const v1(get_camera_pos(), pcenter);
	vector3d v2;
	if (v1.mag() > 1.0E5*min(w1, w2)) return; // too far away
	point pts[4];
	cylinder_quad_projection(pts, cylinder_3dw(p1, p2, w1, w2), v1.get_norm(), v2, npts);
	assert(npts == 4);
	color1.do_glColor();
	
	for (unsigned i = 0; i < 4; ++i) {
		if (i == 2 && color2 != color1) color2.do_glColor();
		glTexCoord2f(((i == 0 || i == 3) ? 0.0 : 1.0), 0.5); // 1D blur texture
		(globalize ? make_pt_global(pts[i]) : pts[i]).do_glVertex();
	}
}


void begin_line_tquad_draw() {

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDisable(GL_LIGHTING);
	enable_blend();
	select_texture(BLUR_TEX);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.01);
	glBegin(GL_QUADS);
}


void end_line_tquad_draw() {

	glEnd();
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
	disable_blend();
	glEnable(GL_LIGHTING);
}


void draw_animated_billboard(point const &pos, float size, float timescale) { // fixed 4x4 animation

	int const frame_id(max(0, min(15, int(16*timescale)))), tx(frame_id&3), ty(frame_id>>2);
	point const camera(get_camera_pos()), gpos(make_pt_global(pos));
	(camera - pos).do_glNormal();
	draw_billboard(gpos, (camera + gpos - pos), up_vector*-1.0, size, size, 0.25*tx, 0.25*ty, 0.25*(tx+1), 0.25*(ty+1)); // upside down
}


void pos_dir_up::draw_frustum() const {

	float const nf_val[2] = {near_, far_};
	point pts[2][4]; // {near, far} x {ll, lr, ur, ul}

	for (unsigned d = 0; d < 2; ++d) {
		point const center(pos + dir*nf_val[d]); // plane center
		float const dy(nf_val[d]*tterm/sqrt(1.0 + A*A)), dx(A*dy); // d*sin(theta)
		pts[d][0] = center - cp*dx - upv_*dy;
		pts[d][1] = center + cp*dx - upv_*dy;
		pts[d][2] = center + cp*dx + upv_*dy;
		pts[d][3] = center - cp*dx + upv_*dy;
	}
	glBegin(GL_QUADS);
	
	// near and far
	for (unsigned d = 0; d < 2; ++d) {
		for (unsigned i = 0; i < 4; ++i) {
			pts[d][i].do_glVertex();
		}
	}

	// sides
	for (unsigned i = 0; i < 4; ++i) {
		pts[0][(i+0)&3].do_glVertex();
		pts[0][(i+1)&3].do_glVertex();
		pts[1][(i+1)&3].do_glVertex();
		pts[1][(i+0)&3].do_glVertex();
	}
	glEnd();
}


int draw_simple_cube(cube_t const &c, bool texture, int in_cur_prim, bool no_normals, int eflags, float texture_scale, vector3d const *const view_dir) {

	if (in_cur_prim != GL_QUADS) {
		if (in_cur_prim >= 0) glEnd();
		glBegin(GL_QUADS);
	}
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if ((eflags & EFLAGS[n][j]) || (view_dir && (((*view_dir)[n] < 0.0) ^ j))) continue; // back facing or disabled
			point pt;
			pt[n] = c.d[n][j];

			if (!no_normals) {
				vector3d norm(zero_vector);
				norm[n] = (2.0*j - 1.0); // -1 or 1
				norm.do_glNormal();
			}
			for (unsigned s = 0; s < 2; ++s) { // d[1] dim
				pt[d[1]] = c.d[d[1]][s];

				for (unsigned k = 0; k < 2; ++k) { // d[0] dim
					pt[d[0]] = c.d[d[0]][k^j^s^1]; // need to orient the vertices differently for each side
						
					if (texture) {
						float const s[2] = {texture_scale*pt[d[1]], texture_scale*pt[d[0]]};
						glTexCoord2fv(s);
					}
					pt.do_glVertex();
				}
			}
		}
	}
	if (in_cur_prim == PRIM_DISABLED) {
		glEnd();
		return in_cur_prim;
	}
	return GL_QUADS;
}


// need to do something with tex coords for scale
void draw_cube(point const &pos, float sx, float sy, float sz, bool texture, unsigned ndiv, bool scale_ndiv,
			   float texture_scale, bool proportional_texture, vector3d const *const view_dir)
{
	point const scale(sx, sy, sz);
	point const xlate(pos - 0.5*scale); // move origin from center to min corner
	float const sizes[3] = {sx, sy, sz};

	if (ndiv == 1) { // non-subdivided version
		glBegin(GL_QUADS);
		
		for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
			unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

			for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
				if (view_dir && (((*view_dir)[n] < 0.0) ^ j)) continue; // back facing
				vector3d norm(zero_vector);
				point pt;
				norm[n] = (2.0*j - 1.0); // -1 or 1
				pt[n]   = j;
				norm.do_glNormal();

				for (unsigned s1 = 0; s1 < 2; ++s1) {
					pt[d[1]] = s1;

					for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
						pt[d[0]] = k^j^s1^1; // need to orient the vertices differently for each side
						
						if (texture) {
							glTexCoord2f((proportional_texture ? sizes[d[1]] : 1.0)*texture_scale*pt[d[1]],
										 (proportional_texture ? sizes[d[0]] : 1.0)*texture_scale*pt[d[0]]);
						}
						(pt*scale + xlate).do_glVertex();
					}
				}
			} // for j
		} // for i
		glEnd();
		return;
	}
	float const step(1.0/float(ndiv));
	unsigned ndivs[3] = {ndiv, ndiv, ndiv};
	float    steps[3] = {step, step, step};

	if (scale_ndiv) {
		float const smax(max(sx, max(sy, sz)));

		for (unsigned i = 0; i < 3; ++i) {
			if (sizes[i] < smax) {
				ndivs[i] = max(1U, unsigned(ceil(ndiv*(sizes[i]/smax))));
				steps[i] = 1.0/float(ndivs[i]);
			}
		}
	}
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (view_dir && (((*view_dir)[n] < 0.0) ^ j)) continue; // back facing
			vector3d norm(zero_vector);
			point pt;
			norm[n] = (2.0*j - 1.0); // -1 or 1
			pt[n]   = j;
			norm.do_glNormal();

			for (unsigned s0 = 0; s0 < ndivs[d[0]]; ++s0) {
				float const va[2] = {steps[d[0]]*(s0 + 1), steps[d[0]]*s0};
				glBegin(GL_QUAD_STRIP);

				for (unsigned s1 = 0; s1 <= ndivs[d[1]]; ++s1) {
					pt[d[1]] = steps[d[1]]*s1;

					for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
						pt[d[0]] = va[k^j]; // need to orient the vertices differently for each side
						
						if (texture) {
							glTexCoord2f((proportional_texture ? sizes[d[1]] : 1.0)*texture_scale*pt[d[1]],
								         (proportional_texture ? sizes[d[0]] : 1.0)*texture_scale*pt[d[0]]);
						}
						(pt*scale + xlate).do_glVertex();
					}
				}
				glEnd();
			}
		} // for j
	} // for i
}


int draw_cylin_quad_proj(cylinder_3dw const &cylin, vector3d const &view_dir, int in_cur_prim, bool no_normals) {

	point pts[4];
	int npts(0);
	vector3d v2; // unused
	cylinder_quad_projection(pts, cylin, view_dir, v2, npts);
	return draw_simple_polygon(pts, npts, view_dir, in_cur_prim, no_normals);
}


int draw_simple_polygon(point const *const points, int npoints, vector3d const &norm, int in_cur_prim, bool no_normals) {

	int prim_type(-1);

	switch (npoints) {
		case 0: return in_cur_prim;
		case 1:  prim_type = GL_POINTS;    break;
		case 2:  prim_type = GL_LINES;     break;
		case 3:  prim_type = GL_TRIANGLES; break;
		case 4:  prim_type = GL_QUADS;     break;
		default: glBegin(GL_POLYGON);
	}
	if (in_cur_prim != prim_type) {
		if (in_cur_prim >= 0) glEnd();
		glBegin(prim_type);
	}
	if (!no_normals) norm.do_glNormal();

	for (int i = 0; i < npoints; ++i) {
		points[i].do_glVertex();
	}
	if (in_cur_prim == PRIM_DISABLED) {
		glEnd();
		return in_cur_prim;
	}
	return prim_type;
}


int draw_simple_extruded_polygon(float thick, point const *const points, int npoints, int in_cur_prim, bool no_normals) {

	assert(points != NULL && (npoints == 3 || npoints == 4));
	thick = fabs(thick);
	vector3d const norm(get_poly_norm(points));
	point pts[2][4];
	gen_poly_planes(points, npoints, norm, thick, pts);
	std::reverse(pts[0], pts[0]+npoints);
	in_cur_prim = draw_simple_polygon(pts[0], npoints, norm*-1, in_cur_prim, no_normals); // draw bottom surface
	std::reverse(pts[0], pts[0]+npoints);
	in_cur_prim = draw_simple_polygon(pts[1], npoints, norm,    in_cur_prim, no_normals); // draw top surface
	
	for (int i = 0; i < npoints; ++i) { // draw sides
		int const ii((i+1)%npoints);
		point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
		in_cur_prim = draw_simple_polygon(side_pts, 4, get_poly_norm(side_pts), in_cur_prim, no_normals);
	}
	return in_cur_prim;
}



void gen_quad_tex_coords(float *tdata, unsigned num, unsigned stride) { // stride is in floats

	for (unsigned i = 0, off = 0; i < num; i++) { // 01 00 10 11 for every quad
		tdata[off] = 0.0; tdata[off+1] = 1.0; off += stride;
		tdata[off] = 0.0; tdata[off+1] = 0.0; off += stride;
		tdata[off] = 1.0; tdata[off+1] = 0.0; off += stride;
		tdata[off] = 1.0; tdata[off+1] = 1.0; off += stride;
	}
}


void gen_quad_tri_tex_coords(float *tdata, unsigned num, unsigned stride) { // stride is in floats

	for (unsigned i = 0, off = 0; i < num; i++) { // for every tri
		tdata[off] = -0.5; tdata[off+1] =  1.0; off += stride;
		tdata[off] =  0.5; tdata[off+1] = -1.0; off += stride;
		tdata[off] =  1.5; tdata[off+1] =  1.0; off += stride;
	}
}


void draw_quads_from_pts(vector<vert_norm> const &points, unsigned draw_num) {

	if (points.empty()) return;
	unsigned const MAX_QUADS(1000);
	unsigned const num(draw_num ? min((unsigned)points.size(), draw_num) : (unsigned)points.size());
	assert((num & 3) == 0);

	if (num+1 < MAX_QUADS) {
		static int init(0);
		static float tdata[8*MAX_QUADS];

		if (!init) {
			gen_quad_tex_coords(tdata, MAX_QUADS, 2);
			init = 1;
		}
		set_array_client_state(1, 1, 1, 0);
		glVertexPointer(3,   GL_FLOAT, sizeof(vert_norm), &points[0].v);
		glNormalPointer(     GL_FLOAT, sizeof(vert_norm), &points[0].n);
		glTexCoordPointer(2, GL_FLOAT, 0, tdata+2); // offset by 2 (one tex coord) to fix texture orientation
		glDrawArrays(GL_QUADS, 0, num);
	}
	else {
		glBegin(GL_QUADS);

		for (unsigned p = 0; p < num; p += 4) { // 00 10 11 01
			for (unsigned i = 0; i < 4; ++i) {
				glTexCoord2f(float(i==1||i==2), float(i==2||i==3));
				points[p+i].n.do_glNormal();
				points[p+i].v.do_glVertex();
			}
		}
		glEnd();
	}
}


// ******************** DLISTS ********************


void free_dlist_block(unsigned *dlists, unsigned num) {

	if (dlists[0] > 0) {
		glDeleteLists(dlists[0], num);
		for (unsigned i = 0; i < num; ++i) dlists[i] = 0;
	}
}


void free_dlists() {

	assert(NUM_RES_DLIST > 0);
	free_dlist_block(predef_dlists, NUM_RES_DLIST);
	free_dlist_block(sphere_dlists, NUM_SPH_DLIST);
	free_dlist_block(cylin_dlists,  NUM_CYLIN_DLIST);
}


inline void draw_half_subdiv_sphere(unsigned ndiv, bool texture) {

	draw_subdiv_sphere_section(all_zeros, 1.0, ndiv, texture, 0.0, 1.0, 0.0, 0.5);
}


void setup_dlists() {

	assert(quadric && NUM_RES_DLIST > 0 && N_SPHERE_DIV > 0);
	
	if (predef_dlists[0] == 0) {
		unsigned const dl0(glGenLists(NUM_RES_DLIST));
		assert(dl0 > 0);

		for (unsigned i = 0; i < NUM_RES_DLIST; ++i) {
			unsigned const dl(dl0 + i);
			assert(glIsList(dl));
			predef_dlists[i] = dl;
			glNewList(dl, GL_COMPILE);
			draw_fast_cylinder(all_zeros, point(0.0, 0.0, 1.0), 1.0, ((i == DLIST_CYLIN || i == DLIST_CYLIN_T) ? 1.0 : 0.0),
				SMALL_NDIV, (i == DLIST_CYLIN_T || i == DLIST_CONE_T));
			glEndList();
		}
	}
	if (sphere_dlists[0] == 0) {
		unsigned const dl0(glGenLists(NUM_SPH_DLIST));
		assert(dl0 > 0);

		for (unsigned i = 1; i <= N_SPHERE_DIV; ++i) {
			for (unsigned tex = 0; tex < 2; ++tex) {
				for (unsigned half = 0; half < 2; ++half) {
					unsigned const index(((i-1) << 2) + (half << 1) + tex), dl(dl0 + index);
					assert(glIsList(dl));
					sphere_dlists[index] = dl;
					glNewList(dl, GL_COMPILE);

					if (half) {
						draw_half_subdiv_sphere(i, (tex != 0));
					}
					else {
						draw_subdiv_sphere(all_zeros, 1.0, i, (tex != 0), 1);
					}
					glEndList();
				}
			}
		}
	}
	if (cylin_dlists[0] == 0) {
		unsigned const dl0(glGenLists(NUM_CYLIN_DLIST));
		assert(dl0 > 0);

		for (unsigned i = 1; i <= NUM_CYLIN_DLIST; ++i) {
			unsigned const dl(dl0 + i-1);
			assert(glIsList(dl));
			cylin_dlists[i-1] = dl;
			glNewList(dl, GL_COMPILE);
			draw_cylinder(1.0, 1.0, 1.0, i, 1, 1);
			glEndList();
		}
	}
}


void draw_cylin_cone_dlist(float r, float l, bool restore_matrix, int type) {

	assert(type < NUM_RES_DLIST);
	if (restore_matrix) glPushMatrix();
	glScalef(r, r, l);
	unsigned const list_id(predef_dlists[type]);
	assert(glIsList(list_id));
	glCallList(list_id);
	if (restore_matrix) glPopMatrix();
}


void draw_cylin_fast(float r1, float r2, float l, int ndiv, bool texture, bool restore_matrix, bool r_approx) {

	if (USE_CYLIN_DLIST && ndiv <= SMALL_NDIV && (r1 == r2 || (r_approx && fabs(r2 - r1) < 0.5*max(r1, r2)))) {
		float const rav(0.5*(r1 + r2)); // inexact
		draw_cylin_cone_dlist(rav, l, restore_matrix, (texture ? DLIST_CYLIN_T : DLIST_CYLIN));
	}
	else if (USE_CYLIN_DLIST && ndiv <= SMALL_NDIV && r2 == 0.0) { // cone pointed up
		draw_cylin_cone_dlist(r1, l, restore_matrix, (texture ? DLIST_CONE_T : DLIST_CONE));
	}
	else {
		draw_fast_cylinder(all_zeros, point(0.0, 0.0, l), r1, r2, ndiv, texture);
	}
}


void draw_sphere_dlist_raw(int ndiv, bool textured, bool half) {

	assert(ndiv <= N_SPHERE_DIV);
	unsigned const dlist_id(sphere_dlists[((ndiv-1) << 2) + (half << 1) + textured]);
	assert(dlist_id > 0);
	//assert(glIsList(dlist_id));
	glCallList(dlist_id);
}


void draw_sphere_dlist(point const &pos, float radius, int ndiv, bool textured, bool half, bool bfc) {

	if (USE_SPHERE_DLIST && ndiv <= N_SPHERE_DIV) { // speedup is highly variable
		assert(ndiv > 0);
		bool const has_xform(radius != 1.0 || pos != all_zeros);
		if (has_xform)        glPushMatrix();
		if (pos != all_zeros) translate_to(pos);
		if (radius != 1.0)    uniform_scale(radius);
		draw_sphere_dlist_raw(ndiv, textured, half);
		if (has_xform)        glPopMatrix();
	}
	else if (half) {
		draw_subdiv_sphere_section(pos, radius, ndiv, textured, 0.0, 1.0, 0.0, 0.5);
	}
	else {
		draw_subdiv_sphere(pos, radius, ndiv, textured, !bfc);
		//draw_sphere_at_tc(pos, radius, ndiv, textured, bfc);
	}
}


void draw_sphere_dlist_back_to_front(point const &pos, float radius, int ndiv, bool textured, bool half) {

	glEnable(GL_CULL_FACE);

	for (unsigned i = 0; i < 2; ++i) { // kind of slow
		glCullFace(i ? GL_BACK : GL_FRONT);
		draw_sphere_dlist(pos, radius, ndiv, textured, half); // cull?, partial sphere?
	}
	glDisable(GL_CULL_FACE);
}


// draw both ends, similar radius, no texture
void draw_rotated_cylinder_dlist(point const &p1, point const &p2, float r, int ndiv, vector3d const &scale) {
	
	assert(ndiv > 0);

	if (!USE_CYLIN_DLIST || ndiv > NUM_CYLIN_DLIST) {
		draw_rotated_cylinder(p1, p2, r, r, ndiv, 1, 1, scale);
		return;
	}
	glPushMatrix();
	xform_cylinder(p1, p2, scale);
	glScalef(r, r, p2p_dist(p1, p2));
	unsigned const list_id(cylin_dlists[ndiv-1]);
	assert(glIsList(list_id));
	glCallList(list_id);
	glPopMatrix();
}



