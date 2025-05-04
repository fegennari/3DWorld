// 3D World - Shape primitive drawing
// by Frank Gennari
// 3/6/06
#include "function_registry.h"
#include "subdiv.h"
#include "upsurface.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "draw_utils.h"
#include "transform_obj.h"


// predefined sphere VBOs
unsigned const MAX_SPHERE_VBO_NDIV = 3*N_SPHERE_DIV/2;
unsigned const NUM_PREDEF_SPHERES  = 4*MAX_SPHERE_VBO_NDIV;

unsigned sphere_vbo_offsets[NUM_PREDEF_SPHERES+1] = {0};
unsigned predef_sphere_vbo(0);
vector_point_norm cylinder_vpn;


extern int display_mode, draw_model;

void setup_shield_shader(shader_t &shader, int noise_tu_id); // ship.cpp

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
	for (unsigned i = 0; i < 2; ++i) {vab[i].normalize();}
}


// create class sd_cylin_d?
// perturb_map usage is untested
vector_point_norm const &gen_cylinder_data(point const ce[2], float radius1, float radius2, unsigned ndiv, vector3d &v12,
										   float const *const perturb_map, float s_beg, float s_end, int force_dim)
{
	assert(ndiv > 0 && ndiv < 100000);
	vector_point_norm &vpn(cylinder_vpn);
	vpn.p.resize(2*ndiv);
	vpn.n.resize(  ndiv);
	float const r[2] = {radius1, radius2};
	// speical case for ndiv==4 as this is a common min LOD size for tree trunk shadows, etc.
	float const css(TWO_PI/(float)ndiv), sin_ds((ndiv == 4) ? 1.0 : sin(css)), cos_ds((ndiv == 4) ? 0.0 : cos(css));
	unsigned s0(0), s1(ndiv);

	if (ce[0].x == ce[1].x && ce[0].y == ce[1].y && s_beg == 0.0 && s_end == 1.0) { // special case optimization for full vertical cylinder
		float const z_sign((ce[1].z > ce[0].z) ? 1.0 : -1.0);
		v12 = z_sign*plus_z;
		float sin_s(0.0), cos_s(1.0);

		for (unsigned S = 0; S < ndiv; ++S) { // build points table
			float const s(sin_s), c(cos_s);
			vpn.p[(S<<1)+0] = ce[0] + vector3d(z_sign*r[0]*s, r[0]*c, 0.0);
			vpn.p[(S<<1)+1] = ce[1] + vector3d(z_sign*r[1]*s, r[1]*c, 0.0);
			sin_s = s*cos_ds + c*sin_ds;
			cos_s = c*cos_ds - s*sin_ds;
		}
	}
	else {
		v12 = (ce[1] - ce[0]).get_norm();
		vector3d vab[2];
		get_ortho_vectors(v12, vab, force_dim);
		s0 = NDIV_SCALE(s_beg); s1 = NDIV_SCALE(s_end);
		if (s1 == ndiv) {s0 = 0;} // make wraparound correct
		s1 = min(ndiv, s1+1); // allow for sn
		float sin_s((s0 == 0) ? 0.0 : sin(s0*css)), cos_s((s0 == 0) ? 1.0 : cos(s0*css));
		// sin(x + y) = sin(x)*cos(y) + cos(x)*sin(y)
		// cos(x + y) = cos(x)*cos(y) - sin(x)*sin(y)

		for (unsigned S = s0; S < s1; ++S) { // build points table
			float const s(sin_s), c(cos_s);
			vpn.p[(S<<1)+0] = ce[0] + (vab[0]*(r[0]*s)) + (vab[1]*(r[0]*c)); // loop unrolled
			vpn.p[(S<<1)+1] = ce[1] + (vab[0]*(r[1]*s)) + (vab[1]*(r[1]*c));
			sin_s = s*cos_ds + c*sin_ds;
			cos_s = c*cos_ds - s*sin_ds;
		}
	}
	if (radius1 == radius2) { // special case optimization for cylinder rather than truncated cone
		float const nscale(1.0/(0.5f*(vpn.p[0] + vpn.p[2]) - ce[0]).mag()); // slightly larger than 1.0/radius
		for (unsigned S = s0; S < s1; ++S) {vpn.n[S] = (0.5f*(vpn.p[S<<1] + vpn.p[((S+1)%ndiv)<<1]) - ce[0])*nscale;} // build normals table
	}
	else {
		bool const npt_r2(radius1 < radius2); // determine normal from longest edge for highest accuracy (required for when r1 == 0.0)
		float nmag_inv(0.0);

		for (unsigned S = s0; S < s1; ++S) { // build normals table
			vector3d const v1(vpn.p[(((S+1)%ndiv)<<1)+npt_r2], vpn.p[(S<<1)+npt_r2]), v2(vpn.p[(S<<1)+1], vpn.p[S<<1]);
			vector3d &normal(vpn.n[S]);
			cross_product(v2, v1, normal);
			if (S == s0) {nmag_inv = 1.0/normal.mag();} // first one (should not have a divide by zero)
			normal *= nmag_inv; //norms[S].normalize();
		}
	}
	if (perturb_map != NULL) { // add in perturbations
		float const ravg(0.5f*(r[0] + r[1]));
		float const pscale[2] = {r[0]/ravg, r[1]/ravg};

		for (unsigned S = s0; S < s1; ++S) {
			for (unsigned i = 0; i < 2; ++i) {vpn.p[(S<<1)+i] += vpn.n[S]*(pscale[i]*perturb_map[S]);}
		}
		// update normals?
	}
	return vpn;
}


// *** sd_sphere_d ***


void sd_sphere_d::gen_points_norms_static(float s_beg, float s_end, float t_beg, float t_end) {

	static sphere_point_norm temp_spn;
	gen_points_norms(temp_spn, s_beg, s_end, t_beg, t_end);
}


void sd_sphere_d::gen_points_norms(sphere_point_norm &cur_spn, float s_beg, float s_end, float t_beg, float t_end) {

	assert(ndiv >= 3 && ndiv <= 512); // sanity check
	bool const is_full(s_beg == 0.0 && s_end == 1.0 && t_beg == 0.0 && t_end == 1.0);

	if (!perturb_map && !surf && cur_spn.points != nullptr && is_full && cur_spn.is_full && ndiv == cur_spn.ndiv) {
		if (radius == cur_spn.radius) { // in this case, the only thing that differs is pos, so we can reuse the normals and simply translate the points
			vector3d const xlate(pos - cur_spn.center);

			for (unsigned s = 0; s < ndiv; ++s) {
				for (unsigned t = 0; t <= ndiv; ++t) {cur_spn.points[s][t] += xlate;}
			}
		}
		else { // in this case, the only thing that differs is pos and radius, so we can reuse the normals and simply translate + scale the points
			float const scale(radius/cur_spn.radius);

			for (unsigned s = 0; s < ndiv; ++s) {
				for (unsigned t = 0; t <= ndiv; ++t) {
					point &p(cur_spn.points[s][t]);
					p -= cur_spn.center; // translate to origin
					p *= scale; // scale the point
					p += pos; // translate to new pos
				}
			}
			cur_spn.radius = radius;
		}
		cur_spn.center = pos;
		points = cur_spn.points;
		norms  = cur_spn.norms;
		return;
	}
	if (cur_spn.points == nullptr || ndiv > cur_spn.ndiv) { // allocate all memory
		if (cur_spn.points != nullptr) {cur_spn.free_data();} // already allocated at least once
		cur_spn.alloc(ndiv);
	}
	else {
		cur_spn.set_pointer_stride(ndiv);
	}
	cur_spn.center  = pos;
	cur_spn.radius  = radius;
	cur_spn.is_full = is_full;
	float const cs_scale(PI/(float)ndiv), cs_scale2(2.0*cs_scale), sin_dt(sin(cs_scale)), cos_dt(cos(cs_scale));
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
			point &pt(cur_spn.points[s][t]);
			float const sv(sin_t), cv(cos_t);
			pt.assign(sv*tvs, sv*tvc, cv); // R*sin(phi)*sin(theta), R*sin(phi)*cos(theta), R*cos(phi)
			sin_t = sv*cos_dt + cv*sin_dt;
			cos_t = cv*cos_dt - sv*sin_dt;
			cur_spn.norms[s][t] = pt; // Note: perturb_map does not affect normals until later
			pt   *= radius;
			pt   += pos;
			if (perturb_map) {pt += cur_spn.norms[s][t]*perturb_map[t+soff];}
			if (surf)        {pt += cur_spn.norms[s][t]*surf->get_height_at(pt);}
		}
	}
	if (perturb_map || surf) { // recalculate vertex/surface normals
		for (unsigned s = s0; s < s1; ++s) {
			for (unsigned t = max(1U, t0); t <= min(t1, ndiv-1); ++t) { // skip the poles
				point const p(cur_spn.points[s][t]);
				vector3d n(zero_vector);

				for (unsigned d = 0; d < 4; ++d) { // s+,t+  t-,s+  s-,t-  t+,s-
					unsigned const si((d==0||d==1) ? (s+1)%ndiv : (s+ndiv-1)%ndiv);
					unsigned const ti((d==0||d==3) ? (t+1)      : (t-1));
					vector3d const nst(cross_product((cur_spn.points[si][t] - p), (cur_spn.points[s][ti] - p)));
					float const nmag(nst.mag());
					if (nmag > TOLERANCE) n += nst*(((d&1) ? -0.25 : 0.25)/nmag);
				}
				cur_spn.norms[s][t] = n;
			}
		}
	}
	points = cur_spn.points;
	norms  = cur_spn.norms;
}


float sd_sphere_d::get_rmax() const { // could calculate this during gen_points_norms

	assert(points);
	float rmax_sq(0.0);

	for (unsigned y = 0; y < ndiv; ++y) {
		for (unsigned x = 0; x <= ndiv; ++x) {
			rmax_sq = max(rmax_sq, p2p_dist_sq(pos, points[y][x]));
		}
	}
	return sqrt(rmax_sq);
}


void sd_sphere_d::set_data(point const &p, float r, int n, float const *pm, float dp, upsurface const *const s) {

	pos         = p;
	radius      = r;
	def_pert    = dp;
	perturb_map = pm;
	surf        = s;
	ndiv        = n;
	assert(radius > 0.0);
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


void sphere_point_norm::free_data() {

	if (points == NULL) {assert(norms == NULL); return;} // already freed (okay)
	assert(points && norms && points[0] && norms[0]);
	matrix_delete_2d(points);
	norms = NULL;
}


// ******************** CYLINDER ********************


void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals, point const &pos, float tscale_s, float tscale_t) {

	assert(r_outer > 0.0);
	bool const disk(r_inner > 0.0);
	vector3d const n(invert_normals ? -plus_z : plus_z);
	float const css((invert_normals ? 1.0 : -1.0)*TWO_PI/(float)ndiv), sin_ds(sin(css)), cos_ds(cos(css));
	float const inner_tscale(r_inner/r_outer);
	float sin_s(0.0), cos_s(1.0);
	static vector<vert_norm_tc> verts;
	if (!disk) {verts.emplace_back(pos, n, 0.5*tscale_s, 0.5*tscale_t);}

	for (unsigned S = 0; S <= (unsigned)ndiv; ++S) {
		float const s(sin_s), c(cos_s);
		if (disk) {verts.emplace_back((pos + point(r_inner*s, r_inner*c, 0.0)), n, 0.5f*(1.0f + inner_tscale*s), (0.5f*(1.0f + inner_tscale*c)));}
		float const ts(0.5*(1.0 + s)*tscale_s), tt(0.5*(1.0 + c)*tscale_t);
		verts.emplace_back((pos + point(r_outer*s, r_outer*c, 0.0)), n, ts, tt);
		sin_s = s*cos_ds + c*sin_ds;
		cos_s = c*cos_ds - s*sin_ds;
	}
	draw_and_clear_verts(verts, (disk ? GL_TRIANGLE_STRIP : GL_TRIANGLE_FAN));
}

void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals, float zval) {
	draw_circle_normal(r_inner, r_outer, ndiv, invert_normals, point(0.0, 0.0, zval));
}


// length can be negative
void draw_cylinder_at(point const &p1, float length, float radius1, float radius2, int ndiv, bool draw_ends, bool first_end_only, bool last_end_only, float tscale_len) {

	assert(ndiv > 0);
	draw_fast_cylinder(p1, p1+vector3d(0.0, 0.0, length), radius1, radius2, ndiv, 1, 0, 0, nullptr, tscale_len); // tex coords?
	if (draw_ends && !last_end_only  && radius1 > 0.0) {draw_circle_normal(0.0, radius1, ndiv, 1, p1);}
	if (draw_ends && !first_end_only && radius2 > 0.0) {draw_circle_normal(0.0, radius2, ndiv, 0, p1+vector3d(0.0, 0.0, length));}
}

void draw_cylinder(float length, float radius1, float radius2, int ndiv, bool draw_ends, bool first_end_only, bool last_end_only, float z_offset, float tscale_len) {
	draw_cylinder_at(point(0.0, 0.0, z_offset), length, radius1, radius2, ndiv, draw_ends, first_end_only, last_end_only, tscale_len);
}


// Note: two_sided_lighting is not entirely correct since it operates on the vertices instead of the faces/fragments
vector3d calc_oriented_normal(point const &p, vector3d const &n, bool two_sided_lighting) {
	return ((two_sided_lighting && dot_product_ptv(n, get_camera_pos(), p) < 0.0) ? -n : n);
}
void create_vert(vert_norm_tc &v, point const &p, vector3d const &n, float ts, float tt, bool two_sided_lighting) {
	v.assign(p, calc_oriented_normal(p, n, two_sided_lighting), ts, tt);
}
void create_vert(vert_norm_texp &v, point const &p, vector3d const &n, texgen_params_t const &tp, bool two_sided_lighting) {
	v = vert_norm_texp(p, calc_oriented_normal(p, n, two_sided_lighting), tp);
}

void gen_cone_triangles(vector<vert_norm_tc> &verts, vector_point_norm const &vpn, bool two_sided_lighting, float tc_t0, float tc_t1, float ts_scale, vector3d const &xlate) {
	
	unsigned const ixoff(verts.size()), ndiv(vpn.n.size());
	verts.resize(3*ndiv + ixoff);
	float const ndiv_inv(1.0/ndiv);

	if (!two_sided_lighting && xlate == zero_vector) { // common case optimization, for example for tree trunks
		for (unsigned s = 0; s < (unsigned)ndiv; ++s) { // Note: always has tex coords
			unsigned const sp((s+ndiv-1)%ndiv), sn((s+1)%ndiv), vix(3*s + ixoff);
			float const ts(ts_scale*(1.0 - s*ndiv_inv));
			verts[vix+0].assign(vpn.p[(s <<1)+1],  vpn.n[s],              ts_scale*0.5,    tc_t1); // one big discontinuity at one position
			verts[vix+1].assign(vpn.p[(sn<<1)+0], (vpn.n[s] + vpn.n[sn]), (ts - ndiv_inv), tc_t0); // normalize?
			verts[vix+2].assign(vpn.p[(s <<1)+0], (vpn.n[s] + vpn.n[sp]), ts,              tc_t0); // normalize?
		}
	}
	else {
		for (unsigned s = 0; s < (unsigned)ndiv; ++s) { // Note: always has tex coords
			unsigned const sp((s+ndiv-1)%ndiv), sn((s+1)%ndiv), vix(3*s + ixoff);
			float const ts(ts_scale*(1.0 - s*ndiv_inv));
			//create_vert(verts[vix+0], vpn.p[(s <<1)+1],        vpn.n[s], (ts - 0.5f*ndiv_inv), tex_scale_len, two_sided_lighting); // small discontinuities at every position
			create_vert(verts[vix+0], vpn.p[(s <<1)+1]+xlate,  vpn.n[s],              ts_scale*0.5,    tc_t1, two_sided_lighting); // one big discontinuity at one position
			create_vert(verts[vix+1], vpn.p[(sn<<1)+0]+xlate, (vpn.n[s] + vpn.n[sn]), (ts - ndiv_inv), tc_t0, two_sided_lighting); // normalize?
			create_vert(verts[vix+2], vpn.p[(s <<1)+0]+xlate, (vpn.n[s] + vpn.n[sp]), ts,              tc_t0, two_sided_lighting); // normalize?
		}
	}
}

void gen_cone_triangles_tp(vector<vert_norm_texp> &verts, vector_point_norm const &vpn, bool two_sided_lighting, texgen_params_t const &tp) {

	unsigned const ixoff(verts.size()), ndiv(vpn.n.size());
	verts.resize(3*ndiv + ixoff);

	for (unsigned s = 0; s < (unsigned)ndiv; ++s) { // Note: always has tex coords
		unsigned const sp((s+ndiv-1)%ndiv), sn((s+1)%ndiv), vix(3*s + ixoff);
		create_vert(verts[vix+0], vpn.p[(s <<1)+1],  vpn.n[s],              tp, two_sided_lighting);
		create_vert(verts[vix+1], vpn.p[(sn<<1)+0], (vpn.n[s] + vpn.n[sn]), tp, two_sided_lighting);
		create_vert(verts[vix+2], vpn.p[(s <<1)+0], (vpn.n[s] + vpn.n[sp]), tp, two_sided_lighting);
	}
}

void gen_cylinder_triangle_strip(vector<vert_norm_tc> &verts, vector_point_norm const &vpn,
	bool two_sided_lighting, float tc_t0, float tc_t1, float ts_scale, vector3d const &xlate)
{
	bool const prev_strip(!verts.empty());
	unsigned const ixoff(prev_strip ? (verts.size() + 2) : 0), ndiv(vpn.n.size()); // 2 extra for connecting with degenerate triangles
	verts.resize(2*(ndiv+1) + ixoff);
	float const ndiv_inv(1.0/ndiv);

	for (unsigned S = 0; S <= ndiv; ++S) { // Note: always has tex coords
		unsigned const s(S%ndiv), vix(2*S + ixoff);
		float const ts(ts_scale*(1.0f - S*ndiv_inv));
		vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
		create_vert(verts[vix+0], vpn.p[(s<<1)+0]+xlate, normal, ts, tc_t0, two_sided_lighting);
		create_vert(verts[vix+1], vpn.p[(s<<1)+1]+xlate, normal, ts, tc_t1, two_sided_lighting);
	}
	if (prev_strip) { // connect previous strip to current strip with degenerate triangles
		verts[ixoff-2] = verts[ixoff-3];
		verts[ixoff-1] = verts[ixoff];
	}
}

void gen_cylinder_quads(vector<vert_norm_tc> &verts, vector_point_norm const &vpn, bool two_sided_lighting, float tex_scale_len) {

	unsigned const ndiv(vpn.n.size());
	float const ndiv_inv(1.0/ndiv);
	unsigned vix(verts.size());
	verts.resize(vix + 4*ndiv);

	for (unsigned i = 0; i < ndiv; ++i) { // Note: always has tex coords
		for (unsigned j = 0; j < 2; ++j) {
			unsigned const S(i + j), s(S%ndiv);
			float const ts(1.0f - S*ndiv_inv);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
			create_vert(verts[vix++], vpn.p[(s<<1)+ j], normal, ts, tex_scale_len*( j), two_sided_lighting);
			create_vert(verts[vix++], vpn.p[(s<<1)+!j], normal, ts, tex_scale_len*(!j), two_sided_lighting);
		}
	} // for i
	assert(vix == verts.size());
}

void gen_cylinder_quads(vector<vert_norm_texp> &verts, vector_point_norm const &vpn, texgen_params_t const &tp, bool two_sided_lighting) {

	unsigned const ndiv(vpn.n.size());

	for (unsigned i = 0; i < ndiv; ++i) { // Note: always has tex coords
		for (unsigned j = 0; j < 2; ++j) {
			unsigned const S(i + j), s(S%ndiv);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
			point const &p1(vpn.p[(s<<1)+!j]), &p2(vpn.p[(s<<1)+ j]);
			verts.emplace_back(p1, calc_oriented_normal(p1, normal, two_sided_lighting), tp);
			verts.emplace_back(p2, calc_oriented_normal(p2, normal, two_sided_lighting), tp);
		}
	} // for i
}


class cylin_vertex_buffer_t {

	bool buffering_enabled;
public:
	vector<vert_norm_tc> tverts, sverts, cverts; // public so they can be accessed from within draw_fast_cylinder()
	vector<vert_wrap_t> verts;

	cylin_vertex_buffer_t() : buffering_enabled(0) {}
	void begin_buffering() {buffering_enabled = 1;}
	void end_cylinder() {if (!buffering_enabled) {draw_and_clear_buffers();}}

	void draw_and_clear_buffers() {
		draw_and_clear_verts(tverts, GL_TRIANGLES);
		draw_and_clear_verts(sverts, GL_TRIANGLE_STRIP);
		buffering_enabled = 0;
	}
};

cylin_vertex_buffer_t cylin_vertex_buffer;

void begin_cylin_vertex_buffering() {cylin_vertex_buffer.begin_buffering();}
void flush_cylin_vertex_buffer   () {cylin_vertex_buffer.draw_and_clear_buffers();}

// draw_sides_ends: 0 = draw sides only, 1 = draw sides and ends, 2 = draw ends only, 3 = pt1 end, 4 = pt2 end
void add_cylin_ends(float radius1, float radius2, int ndiv, bool texture, int draw_sides_ends,
	vector3d const &v12, point const ce[2], point const &xlate, vector_point_norm const &vpn)
{
	if (draw_sides_ends == 0) return;
	cylin_vertex_buffer_t &cvb(cylin_vertex_buffer);
	float const theta_mult(TWO_PI/ndiv);

	for (unsigned i = 0; i < 2; ++i) {
		if ((i ? radius2 : radius1) == 0.0 || (draw_sides_ends == 3+(!i))) continue;
		vector3d const normal(i ? v12 : -v12);
		cvb.cverts.resize(ndiv+2);
		cvb.cverts[0].assign(ce[i]+xlate, normal, 0.5, 0.5);

		for (unsigned S = 0; S <= (unsigned)ndiv; ++S) {
			unsigned const ss(S%ndiv), s(i ? (ndiv - ss - 1) : ss);
			float ts(0.0), tt(0.0);

			if (texture) { // inefficient, but uncommon
				float const theta(theta_mult*s);
				ts = 0.5*(1.0 + sinf(theta));
				tt = 0.5*(1.0 + cosf(theta));
			}
			cvb.cverts[S+1].assign(vpn.p[(s<<1)+i]+xlate, normal, ts, tt);
		}
		draw_and_clear_verts(cvb.cverts, GL_TRIANGLE_FAN); // triangle fans can't be buffered
	}
}

// draw_sides_ends: 0 = draw sides only, 1 = draw sides and ends, 2 = draw ends only, 3 = sides + pt1 end, 4 = sides + pt2 end
void draw_fast_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, bool texture, int draw_sides_ends,
	bool two_sided_lighting, float const *const perturb_map, float tex_scale_len, float tex_t_start, point const *inst_pos, unsigned num_insts, float tex_width_scale)
{
	assert(radius1 > 0.0 || radius2 > 0.0);
	point const ce[2] = {p1, p2};
	vector3d v12; // (ce[1] - ce[0]).get_norm()
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, ndiv, v12, perturb_map));
	cylin_vertex_buffer_t &cvb(cylin_vertex_buffer);
	if (inst_pos == NULL) {inst_pos = &zero_vector; num_insts = 1;} // default identity transform
	assert(num_insts > 0);

	for (unsigned inst = 0; inst < num_insts; ++inst) {
		if (draw_sides_ends == 2) {
			// draw ends only - nothing to do here
		}
		else if (radius2 == 0.0) { // cone (Note: still not perfect for pine tree trunks and enforcer ships)
			gen_cone_triangles(cvb.tverts, vpn, two_sided_lighting, tex_t_start, tex_scale_len+tex_t_start, tex_width_scale, inst_pos[inst]); // triangles
		}
		else {
			gen_cylinder_triangle_strip(cvb.sverts, vpn, two_sided_lighting, tex_t_start, tex_scale_len+tex_t_start, tex_width_scale, inst_pos[inst]); // triangle strip
		}
		if (draw_sides_ends != 0) {add_cylin_ends(radius1, radius2, ndiv, texture, draw_sides_ends, v12, ce, inst_pos[inst], vpn);} // triangle fan; Note: TSL doesn't apply here
	} // for inst
	cvb.end_cylinder();
}

// no normals or tex coords
void draw_shadow_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, int draw_ends, float const *const perturb_map) {

	assert(radius1 > 0.0 || radius2 > 0.0);
	point const ce[2] = {p1, p2};
	vector3d v12; // (ce[1] - ce[0]).get_norm()
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, ndiv, v12, perturb_map));
	vector<vert_wrap_t> &verts(cylin_vertex_buffer.verts);
	verts.resize(2*(ndiv+1));

	for (unsigned S = 0; S <= (unsigned)ndiv; ++S) {
		unsigned const s(S%ndiv), vix(2*S);
		verts[vix+0].v = vpn.p[(s<<1)+0];
		verts[vix+1].v = vpn.p[(s<<1)+1];
	}
	draw_and_clear_verts(verts, GL_TRIANGLE_STRIP);
	if (!draw_ends) return;

	for (unsigned i = 0; i < 2; ++i) {
		if ((i ? radius2 : radius1) == 0.0) continue;
		verts.resize(ndiv+2);
		verts[0].v = ce[i];

		for (unsigned S = 0; S <= (unsigned)ndiv; ++S) {
			unsigned const ss(S%ndiv), s(i ? (ndiv - ss - 1) : ss);
			verts[S+1].v = vpn.p[(s<<1)+i];
		}
		draw_and_clear_verts(verts, GL_TRIANGLE_FAN); // triangle fans can't be buffered
	}
}


void draw_cylin_fast(float r1, float r2, float l, int ndiv, bool texture, float tex_scale_len, float z_offset) {
	draw_fast_cylinder(point(0.0, 0.0, z_offset), point(0.0, 0.0, z_offset+l), r1, r2, ndiv, texture, 0, 0, NULL, tex_scale_len);
}


void draw_cylindrical_section(float length, float r_inner, float r_outer, int ndiv, bool texture, float tex_scale_len, float z_offset) {

	assert(r_outer > 0.0 && r_inner >= 0.0 && length >= 0.0 && ndiv > 0 && r_outer >= r_inner);
	draw_cylin_fast(r_outer, r_outer, length, ndiv, texture, tex_scale_len, z_offset);

	if (r_inner != r_outer) {
		if (r_inner != 0.0) {draw_cylin_fast(r_inner, r_inner, length, ndiv, texture, tex_scale_len, z_offset);}
		draw_circle_normal(r_inner, r_outer, ndiv, 1, z_offset);
		draw_circle_normal(r_inner, r_outer, ndiv, 0, z_offset+length);
	}
}


// ******************** SPHERE ********************


struct sphere_verts_t {
	vector<vert_norm> vn;
	vector<vert_norm_tc> vntc;

	void draw(bool texture, bool use_quads=0, bool do_clear=1) {
		if (use_quads) {
			if (texture) {draw_quad_verts_as_tris(vntc);} else {draw_quad_verts_as_tris(vn);}
		}
		else {
			if (texture) {draw_verts(vntc, GL_TRIANGLE_STRIP);} else {draw_verts(vn, GL_TRIANGLE_STRIP);}
		}
		if (do_clear) {vn.clear(); vntc.clear();}
	}
};

// back face culling requires the sphere to be untransformed (or the vertices to be per-transformed)
void sd_sphere_d::draw_subdiv_sphere(point const &vfrom, int texture, bool disable_bfc, unsigned char const *const render_map,
									 float const *const exp_map, point const *const pt_shift, float expand,
									 float s_beg, float s_end, float t_beg, float t_end, bool back_to_front) const
{
	assert(!render_map || disable_bfc);
	assert(ndiv > 0);
	float const ndiv_inv(1.0/float(ndiv)), rv(def_pert + radius), rv_sq(rv*rv), tscale(texture);
	float const toler(1.0E-6*radius*radius + rv_sq), dmax(rv + 0.1*radius), dmax_sq(dmax*dmax);
	bool const use_quads(render_map || pt_shift || exp_map || expand != 0.0);
	if (expand != 0.0) {expand *= 0.25;} // 1/4, for normalization
	unsigned const s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end)), t0(NDIV_SCALE(t_beg)), t1(NDIV_SCALE(t_end));
	static sphere_verts_t sv;

	for (unsigned s = s0; s < s1; ++s) {
		unsigned const sn((s+1)%ndiv), snt(s+1);

		if (use_quads) { // use slower quads - no back face culling
			for (unsigned t = t0; t < t1; ++t) {
				unsigned const ix(s*(ndiv+1)+t), tn(t+1);
				if (render_map && !render_map[ix]) continue;
				float const exp(expand + (exp_map ? exp_map[ix] : 0.0));
				point          pts[4]     = {points[s][t], points[sn][t], points[sn][tn], points[s][tn]};
				vector3d const normals[4] = {norms [s][t], norms [sn][t], norms [sn][tn], norms [s][tn]};
				
				if (exp != 0.0) { // average the normals
					vector3d const quad_norm(normals[0] + normals[1] + normals[2] + normals[3]);
					for (unsigned i = 0; i < 4; ++i) {pts[i] += quad_norm*exp;}
				}
				for (unsigned i = 0; i < 4; ++i) {
					if (pt_shift) {pts[i] += pt_shift[ix];}
					if (texture) {sv.vntc.emplace_back(pts[i], normals[i], tscale*(1.0f - (((i&1)^(i>>1)) ? snt : s)*ndiv_inv), tscale*(1.0f - ((i>>1) ? tn : t)*ndiv_inv));}
					else {sv.vn.emplace_back(pts[i], normals[i]);}
				}
			} // for t
		}
		else { // use triangle strip
			if (s != s0) { // add degenerate triangles to preserve the triangle strip (only slightly faster than using multiple triangle strips)
				for (unsigned d = 0; d < 2; ++d) {
					unsigned const T(d ? t0 : t1);
					if (texture) {sv.vntc.emplace_back(points[s][T], norms[s][T], 0.0, 0.0);}
					else {sv.vn.emplace_back(points[s][T], norms[s][T]);}
				}
			}
			float const tc1(tscale*(1.0f - s*ndiv_inv)), tc2(tscale*(1.0f - snt*ndiv_inv));

			for (unsigned t = t0; t <= t1; ++t) {
				if (!disable_bfc) {
					bool draw(0);

					for (unsigned d = 0; d < 2 && !draw; ++d) {
						point    const &pt(d ? points[sn][t] : points[s][t]);
						vector3d const &n (d ? norms [sn][t] : norms [s][t]);
						float const dp(dot_product_ptv(n, vfrom, pt));
						if (dp >= 0.0) {draw = 1;}
						else if (perturb_map != NULL) { // sort of a hack, not always correct but good enough in most cases
							float const dist_sq(p2p_dist_sq(pos, pt));
							if (dist_sq > toler && (dist_sq > dmax_sq || dp > -0.3*p2p_dist(vfrom, pt))) {draw = 1;}
						}
					}
					if (!draw) {continue;}
				}
				if (!texture) {
					sv.vn.emplace_back(points[s ][t], norms[s ][t]);
					sv.vn.emplace_back(points[sn][t], norms[sn][t]);
				}
				else {
					float const ts(tscale*(1.0f - t*ndiv_inv));
					sv.vntc.emplace_back(points[s ][t], norms[s ][t], tc1, ts);
					sv.vntc.emplace_back(points[sn][t], norms[sn][t], tc2, ts);
				}
			} // for t
		}
	} // for s
	if (back_to_front) {
		glCullFace(GL_FRONT);
		sv.draw(texture, use_quads, 0); // do_clear=0
		glCullFace(GL_BACK);
	}
	sv.draw(texture, use_quads);
}


void draw_cube_mapped_sphere(point const &center, float radius, unsigned ndiv, bool texture) {

	assert(radius > 0.0 && ndiv > 0);
	float const tstep(1.0/ndiv), vstep(2.0*tstep);
	static sphere_verts_t sv;

	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d1(i), d2((i+1)%3), dn((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			point pt;
			pt[dn] = (j ? 1.0 : -1.0);

			for (unsigned s = 0; s < ndiv; ++s) {
				pt[d1] = -1.0f + s*vstep;

				for (unsigned T = 0; T <= ndiv; ++T) {
					unsigned const t(j ? T : ndiv-T); // reverse between sides
					pt[d2] = -1.0f + t*vstep;
					point pt2(pt);

					for (unsigned k = 0; k < 2; ++k) {
						vector3d const n(pt2.get_norm());
						point const pos(center + n*radius);
						if (texture) {sv.vntc.emplace_back(pos, n, (s+k)*tstep, t*tstep);} else {sv.vn.emplace_back(pos, n);}
						pt2[d1] += vstep;
					}
				} // for t
				sv.draw(texture);
			} // for s
		} // for j
	} // for i
}


// used for scenery, not using vertex_type_t here
void sd_sphere_d::get_quad_points(vector<vert_norm_tc> &quad_pts, vector<unsigned> *indices, bool use_tri_strip, float s_beg, float s_end, float t_beg, float t_end) const {

	assert(ndiv > 0);
	float const ndiv_inv(1.0/float(ndiv));
	bool const is_full(s_beg == 0.0 && s_end == 1.0 && t_beg == 0.0 && t_end == 1.0);
	if (is_full && !use_tri_strip && indices == nullptr && quad_pts.empty()) {quad_pts.reserve(4*ndiv*ndiv);}
	unsigned const s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end)), t0(NDIV_SCALE(t_beg)), t1(NDIV_SCALE(t_end)), stride(t1 - t0 + 1);
	unsigned const s_stop(s1 + (indices != nullptr && !use_tri_strip)); // one extra s iteration for indexed quads to get the final texture coord
	
	for (unsigned s = s0; s < s_stop; ++s) {
		unsigned const sn((s+1)%ndiv), snt(min((s+1), ndiv));

		if (use_tri_strip && indices != nullptr) { // indexed triangle strip mode
			// add degenerate triangle to preserve the triangle strip, even for the first strip in case another sphere is merged into the same vert stream (leafy plants)
			indices->push_back(quad_pts.size());
			indices->push_back(quad_pts.size()+stride-1);

			for (unsigned t = t0; t <= t1; ++t) {
				unsigned const ix(quad_pts.size());
				quad_pts.emplace_back(points[s][t], norms[s][t], (1.0f - s*ndiv_inv), (1.0f - t*ndiv_inv));
				indices->push_back(ix);
				indices->push_back(ix+stride); // +s
			}
		}
		else if (use_tri_strip) { // triangle strip mode
			// add degenerate triangle to preserve the triangle strip, even for the first strip in case another sphere is merged into the same vert stream (leafy plants)
			for (unsigned d = 0; d < 2; ++d) {quad_pts.emplace_back(points[s][d ? t0 : t1], norms[s][d ? t0 : t1], 0, 0);}

			for (unsigned t = t0; t <= t1; ++t) {
				quad_pts.emplace_back(points[s ][t], norms[s ][t], (1.0f - s  *ndiv_inv), (1.0f - t*ndiv_inv));
				quad_pts.emplace_back(points[sn][t], norms[sn][t], (1.0f - snt*ndiv_inv), (1.0f - t*ndiv_inv));
			}
		}
		else if (indices != nullptr) { // indexed quads/triangles mode
			for (unsigned t = t0; t <= t1; ++t) {
				unsigned const ix(quad_pts.size()), S(s%ndiv);
				quad_pts.emplace_back(points[S][t], norms[S][t], (1.0f - s*ndiv_inv), (1.0f - t*ndiv_inv));
				if (t == t1 || s == s1) continue; // no indices added for last s or t values
				indices->push_back(ix);
				indices->push_back(ix+stride); // +s
				indices->push_back(ix+stride+1); // +s +t
				indices->push_back(ix+1); // +t
			} // for t
		}
		else { // quads mode
			for (unsigned t = t0; t < t1; ++t) {
				point          pts[4]     = {points[s][t], points[sn][t], points[sn][t+1], points[s][t+1]};
				vector3d const normals[4] = {norms [s][t], norms [sn][t], norms [sn][t+1], norms [s][t+1]};

				for (unsigned i = 0; i < 4; ++i) {
					quad_pts.emplace_back(pts[i], normals[i], (1.0f - (((i&1)^(i>>1)) ? snt : s)*ndiv_inv), (1.0f - ((i>>1) ? t+1 : t)*ndiv_inv));
				}
			} // for t
		}
	} // for s
}


void sd_sphere_d::get_triangles(vector<vert_wrap_t> &verts) const {

	assert(ndiv > 0);
	
	for (unsigned s = 0; s < ndiv; ++s) {
		unsigned const sn((s+1)%ndiv);

		for (unsigned t = 0; t < ndiv; ++t) {
			verts.push_back(points[s ][t  ]); // 0
			verts.push_back(points[sn][t  ]); // 1
			verts.push_back(points[sn][t+1]); // 2
			verts.push_back(points[s ][t  ]); // 0
			verts.push_back(points[sn][t+1]); // 2
			verts.push_back(points[s ][t+1]); // 3
		}
	}
}


void sd_sphere_d::get_triangle_strip_pow2(vector<vertex_type_t> &verts, unsigned skip, float s_beg, float s_end, float t_beg, float t_end) const {

	assert(ndiv > 0);
	float const ndiv_inv(1.0/float(ndiv));
	unsigned const s0(NDIV_SCALE(s_beg)), s1(NDIV_SCALE(s_end)), t0(NDIV_SCALE(t_beg)), t1(NDIV_SCALE(t_end));

	for (unsigned s = s0; s < s1; s += skip) {
		s = min(s, ndiv-1);
		unsigned const sn((s+skip)%ndiv), snt(min((s+skip), ndiv));

		if (s != s0) { // add degenerate triangle to preserve the triangle strip
			for (unsigned d = 0; d < 2; ++d) {verts.emplace_back(points[s][d ? t0 : t1], norms[s][d ? t0 : t1], 0, 0);}
		}
		for (unsigned t = t0; t <= t1; t += skip) {
			t = min(t, ndiv);
			verts.emplace_back(points[s ][t], norms[s ][t], (1.0f - s  *ndiv_inv), (1.0f - t*ndiv_inv));
			verts.emplace_back(points[sn][t], norms[sn][t], (1.0f - snt*ndiv_inv), (1.0f - t*ndiv_inv));
		}
	} // for s
}


void sd_sphere_d::get_triangle_vertex_list(vector<vertex_type_t> &verts) const {

	float const ndiv_inv(1.0/float(ndiv));

	for (unsigned s = 0; s <= ndiv; ++s) {
		unsigned const six(s%ndiv);

		for (unsigned t = 0; t <= ndiv; ++t) {
			verts.emplace_back(points[six][t], norms[six][t], (1.0f - s*ndiv_inv), (1.0f - t*ndiv_inv));
		}
	}
	//assert(verts.size() < (1ULL << 8*sizeof(index_type_t)));
}


void sd_sphere_d::get_triangle_index_list_pow2(vector<index_type_t> &indices, unsigned skip) const {

	unsigned const stride(ndiv + 1);

	for (unsigned s = 0; s < ndiv; s += skip) {
		if (s != 0) { // add degenerate triangle to preserve the triangle strip
			for (unsigned d = 0; d < 2; ++d) {indices.push_back(s*stride);}
		}
		for (unsigned t = 0; t <= ndiv; t += skip) {
			t = min(t, ndiv);
			indices.push_back((s+0)   *stride + t);
			indices.push_back((s+skip)*stride + t);
		}
	}
}


void sd_sphere_d::get_faceted_triangles(vector<vertex_type_t> &verts) const {

	assert(ndiv > 0);
	float const ndiv_inv(1.0/float(ndiv));
	
	for (unsigned s = 0; s < ndiv; ++s) {
		unsigned const sn((s+1)%ndiv);

		for (unsigned t = 0; t < ndiv; ++t) {
			unsigned const sixs[2][3] = {{s, sn, sn}, {s, sn, s}}, tixs[2][3] = {{t, t, t+1}, {t, t+t, t+1}};
			triangle const tris[2] = {triangle(points[s][t], points[sn][t  ], points[sn][t+1]),
				                      triangle(points[s][t], points[sn][t+1], points[s ][t+1])};
			for (unsigned d = 0; d < 2; ++d) {
				vector3d const normal(tris[d].get_normal()); // face normal
				UNROLL_3X(verts.emplace_back(tris[d].pts[i_], normal, (1.0f - sixs[d][i_]*ndiv_inv), (1.0f - tixs[d][i_]*ndiv_inv)););
			}
		}
	}
}


void sd_sphere_vbo_d::ensure_vbos() {

	// Note: have to re-bind VBO and call vertex_type_t::set_vbo_arrays() during draw_setup() each time because it's drawn with different shaders (asteroids and comets)
	if (!vbo) {
		assert(!ivbo);
		vector<vertex_type_t> verts;
		vector<index_type_t> indices;
		if (faceted) {get_faceted_triangles(verts);} else {get_triangle_vertex_list(verts);}

		if (!faceted) {
			assert(ix_offsets.empty());
			ix_offsets.push_back(0);

			for (unsigned n = ndiv, skip = 1, ix_ix = 0; n >= 4; n >>= 1, skip <<= 1, ++ix_ix) {
				get_triangle_index_list_pow2(indices, skip);
				ix_offsets.push_back(indices.size());
			}
		}
		create_and_upload(verts, indices);
	}
	assert(faceted || !ix_offsets.empty());
}

void sd_sphere_vbo_d::clear_vbos() {
	indexed_vao_manager_t::clear_vbos();
	ix_offsets.clear();
}

unsigned calc_lod_pow2(unsigned max_ndiv, unsigned ndiv) {

	ndiv = max(ndiv, 4U);
	unsigned lod(0);
	for (unsigned n = (max_ndiv >> 1); ndiv <= n; n >>= 1, ++lod) {}
	return lod;
}

unsigned sd_sphere_vbo_d::draw_setup(unsigned draw_ndiv) {

	ensure_vbos();
	pre_render(!faceted, 1); // do_bind_vbo=1
	vertex_type_t::set_vbo_arrays();
	return (faceted ? 0 : calc_lod_pow2(ndiv, draw_ndiv));
}

void sd_sphere_vbo_d::draw_ndiv_pow2_vbo(unsigned draw_ndiv) {

	unsigned const lod(draw_setup(draw_ndiv));
	if (faceted) {glDrawArrays(GL_TRIANGLES, 0, 6*ndiv*ndiv);} // ndiv is ignored
	else {glDrawRangeElements(GL_TRIANGLE_STRIP, 0, (ndiv+1)*(ndiv+1), get_count(lod), get_index_type_enum(), get_index_ptr(lod));}
	++num_frame_draw_calls;
	post_render();
}

void sd_sphere_vbo_d::draw_instances(unsigned draw_ndiv, instance_render_t &inst_render) {

	unsigned const lod(draw_setup(draw_ndiv));
	if (faceted) {inst_render.draw_and_clear(GL_TRIANGLES, 6*ndiv*ndiv, vbo);} // ndiv is ignored
	else {inst_render.draw_and_clear(GL_TRIANGLE_STRIP, get_count(lod), vbo, get_index_type_enum(), get_index_ptr(lod));}
	post_render();
}


void get_sphere_triangles(vector<vert_wrap_t> &verts, point const &pos, float radius, int ndiv) {

	sd_sphere_d sd(pos, radius, ndiv);
	sd.gen_points_norms_static();
	sd.get_triangles(verts);
}

void add_sphere_quads(vector<vert_norm_tc> &verts, vector<unsigned> *indices, point const &pos, float radius, int ndiv,
	bool use_tri_strip, float s_beg, float s_end, float t_beg, float t_end)
{
	sd_sphere_d sd(pos, radius, ndiv);
	sd.gen_points_norms_static(s_beg, s_end, t_beg, t_end);
	sd.get_quad_points(verts, indices, use_tri_strip, s_beg, s_end, t_beg, t_end);
}


// non-collision object version
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, point const &vfrom, float const *perturb_map,
						int texture, bool disable_bfc, unsigned char const *const render_map, float const *const exp_map,
						point const *const pt_shift, float expand, float s_beg, float s_end, float t_beg, float t_end)
{
	sd_sphere_d sd(pos, radius, ndiv, perturb_map);
	sd.gen_points_norms_static(s_beg, s_end, t_beg, t_end);
	sd.draw_subdiv_sphere(vfrom, texture, disable_bfc, render_map, exp_map, pt_shift, expand, s_beg, s_end, t_beg, t_end);
}

void draw_subdiv_sphere(point const &pos, float radius, int ndiv, int texture, bool disable_bfc) {
	draw_subdiv_sphere(pos, radius, ndiv, (disable_bfc ? all_zeros : get_camera_all()), NULL, texture, disable_bfc);
}

void draw_subdiv_sphere_section(point const &pos, float radius, int ndiv, int texture, float s_beg, float s_end, float t_beg, float t_end) {
	draw_subdiv_sphere(pos, radius, ndiv, all_zeros, NULL, texture, 1, NULL, NULL, NULL, 0.0, s_beg, s_end, t_beg, t_end);
}


void rotate_sphere_tex_to_dir(vector3d const &dir) { // dir must be normalized
	fgRotate(atan2(-dir.y, -dir.x)*TO_DEG-90.0, 0.0, 0.0, 1.0); // negate because dir was backwards
	fgRotate(asinf(-dir.z)*TO_DEG, 1.0, 0.0, 0.0);
}


void draw_single_colored_sphere(point const &pos, float radius, int ndiv, colorRGBA const &color) {

	shader_t s;
	s.begin_color_only_shader(color);
	draw_sphere_vbo(pos, radius, ndiv, 0);
	s.end_shader();
}


// ******************** TORUS ********************


vector<float> const &gen_torus_sin_cos_vals(unsigned ndivi) {

	assert(ndivi > 0);
	static vector<float> sin_cos;
	if (sin_cos.size() == 2*ndivi) return sin_cos; // since sin_cos only depends on ndivi, we only need to recompute it when ndivi changes
	sin_cos.resize(2*ndivi);
	float const dt(TWO_PI/ndivi);

	for (unsigned t = 0; t < ndivi; ++t) {
		float const phi(t*dt);
		sin_cos[(t<<1)+0] = cos(phi);
		sin_cos[(t<<1)+1] = sin(phi);
	}
	return sin_cos;
}

// always textured
void draw_rot_torus(point const &center, vector3d const &dir, float ri, float ro, unsigned ndivi, unsigned ndivo, float tex_scale_i, float tex_scale_o) {

	assert(ndivi > 2 && ndivo > 2);
	float const ts(tex_scale_o/ndivo), tt(tex_scale_i/ndivi), ds(TWO_PI/ndivo), cds(cos(ds)), sds(sin(ds));
	static vector<vert_norm_tc> verts;
	verts.resize(2*(ndivi+1), vert_norm_tc(all_zeros, zero_vector, 0.0, 0.0));
	vector<float> const &sin_cos(gen_torus_sin_cos_vals(ndivi));
	vector3d vab[2];
	get_ortho_vectors(dir, vab);
	
	for (unsigned s = 0; s < ndivo; ++s) { // outer
		float const theta(s*ds), ct(cos(theta)), st(sin(theta)), ct2(ct*cds - st*sds), st2(st*cds + ct*sds);
		point const pos [2] = {(vab[0]*ct + vab[1]*st), (vab[0]*ct2 + vab[1]*st2)};
		point const vpos[2] = {(center + pos[0]*ro), (center + pos[1]*ro)};

		for (unsigned t = 0; t <= ndivi; ++t) { // inner
			unsigned const t_((t == ndivi) ? 0 : t);
			float const cp(sin_cos[(t_<<1)+0]), sp(sin_cos[(t_<<1)+1]);

			for (unsigned i = 0; i < 2; ++i) {
				vector3d const delta(pos[1-i]*sp + dir*cp);
				verts[(t<<1)+i].assign((vpos[1-i] + delta*ri), delta, ts*(s+1-i), tt*t);
			}
		} // for t
		draw_verts(verts, GL_TRIANGLE_STRIP);
	} // for s
}

// in z-plane, always textured
void draw_torus(point const &center, float ri, float ro, unsigned ndivi, unsigned ndivo, float tex_scale_i, float tex_scale_o) {
	draw_rot_torus(center, plus_z, ri, ro, ndivi, ndivo, tex_scale_i, tex_scale_o);
}


// ******************** QUAD/CUBE/POLYGON ********************


void rotate_towards_camera(point const &pos) {
	rotate_into_plus_z((get_camera_pos() - pos));
}


void enable_flares(int tid) { // used for clouds and smoke

	glDepthMask(GL_FALSE); // not quite right - prevents flares from interfering with each other but causes later shapes to be drawn on top of the flares
	enable_blend();
	if (draw_model == 0) {select_texture((tid < 0) ? BLUR_TEX : tid);}
}

void disable_flares() {

	disable_blend();
	glDepthMask(GL_TRUE);
}


void draw_one_tquad(float x1, float y1, float x2, float y2, float z, int prim_type) { // Note: normal is +z

	vert_norm_tc verts[4];
	verts[0] = vert_norm_tc(point(x1, y1, z), plus_z, 0, 0); // clockwise
	verts[1] = vert_norm_tc(point(x1, y2, z), plus_z, 0, 1);
	verts[2] = vert_norm_tc(point(x2, y2, z), plus_z, 1, 1);
	verts[3] = vert_norm_tc(point(x2, y1, z), plus_z, 1, 0);
	draw_verts(verts, 4, prim_type); // GL_TRIANGLE_FAN (quads) or GL_PATCHES
}

void draw_tquad(float xsize, float ysize, float z, int prim_type) { // Note: normal is +z
	draw_one_tquad(-xsize, -ysize, xsize, ysize, z, prim_type);
}


// ordered p1+, p1-, p2-, p2+
int get_line_as_quad_pts(point const &p1, point const &p2, float w1, float w2, point pts[4]) {

	int npts(0);
	vector3d const v1(get_camera_pos(), (p1 + p2)*0.5);
	float const dmax(1.0E5*max(w1, w2));
	if (v1.mag_sq() > dmax*dmax) return 0; // too far away
	cylinder_quad_projection(pts, p1, p2, w1, w2, v1, npts);
	assert(npts == 3 || npts == 4);
	return npts;
}


void line_tquad_draw_t::add_line_as_tris(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2,
	point const* const prev, point const *const next, bool make_global, bool use_2d_tex_coords)
{
	//assert(color1.is_valid() && color2.is_valid()); // validate
	if (prev || next) {assert(w1 > 0.0 && w2 > 0.0);} else {assert(w1 >= 0.0 && w2 >= 0.0);}
	point pts[5];
	int const npts(get_line_as_quad_pts(p1, p2, w1, w2, pts));
	float const tc1(use_2d_tex_coords ? 0.0 : 0.5), tc2(use_2d_tex_coords ? 1.0 : 0.5);
	if (npts == 0) return;

	if (npts == 3) { // single triangle
		assert(!prev && !next);
		if (make_global) {UNROLL_3X(pts[i_] = make_pt_global(pts[i_]);)}
		verts.emplace_back(pts[0], 0.0, tc1, color1);
		verts.emplace_back(pts[1], ((w1 == 0.0) ? 1.0 : 0.0), ((w1 == 0.0) ? tc2 : tc1), ((w1 == 0.0) ? color2 : color1));
		verts.emplace_back(pts[2], 1.0, tc2, color2);
		return;
	}
	if (prev && *prev != p1) {
		vector3d const dv(w1*cross_product((get_camera_pos() - (*prev + p1)*0.5), (p1 - *prev)).get_norm());
		pts[0] = 0.5*(pts[0] + (p1 + dv)); // average the points
		pts[1] = 0.5*(pts[1] + (p1 - dv)); // average the points
	}
	if (next && *next != p2) {
		vector3d const dv(w2*cross_product((get_camera_pos() - (p2 + *next)*0.5), (*next - p2)).get_norm());
		pts[2] = 0.5*(pts[2] + (p2 - dv)); // average the points
		pts[3] = 0.5*(pts[3] + (p2 + dv)); // average the points
	}
	pts[4] = p2;
	if (make_global) {for (unsigned i = 0; i < 5; ++i) {pts[i] = make_pt_global(pts[i]);}}
	color_wrapper cw1, cw2; cw1.set_c4(color1); cw2.set_c4(color2);
	verts.emplace_back(pts[2], 0.0, tc2, cw2.c);
	verts.emplace_back(pts[1], 0.0, tc1, cw1.c);
	verts.emplace_back(pts[4], 0.5, tc2, cw2.c);
	verts.push_back(verts.back()); // duplicate
	verts.emplace_back(pts[1], 0.0, tc1, cw1.c);
	verts.emplace_back(pts[0], 1.0, tc1, cw1.c);
	verts.push_back(verts.back()); // duplicate
	verts.emplace_back(pts[3], 1.0, tc2, cw2.c);
	verts.emplace_back(pts[4], 0.5, tc2, cw2.c);
}


void line_tquad_draw_t::draw(float noise_scale) const { // supports quads and triangles

	if (empty()) return;
	shader_t s;

	if (0) { // probably more efficient, and preserves sharp/bright lines
		s.set_vert_shader("no_lighting_tex_coord");
		s.set_frag_shader("line_draw_halo");
		s.begin_shader();
	}
	else if (noise_scale > 0.0) {
		s.set_prefix("#define LINE_MODE", 1); // FS
		setup_shield_shader(s, 1);
		s.add_uniform_float("noise_scale", noise_scale);
		select_texture((draw_model != 0) ? (int)WHITE_TEX : (int)BLUR_TEX);
	}
	else { // texture mipmaps perform antialiasing on distant lines, which looks nice
		s.begin_simple_textured_shader(0.01);
		select_texture((draw_model != 0) ? (int)WHITE_TEX : (int)BLUR_TEX);
	}
	enable_blend();
	draw_tri_verts();
	disable_blend();
	s.end_shader();
}


void pos_dir_up::draw_frustum() const {

	point pts[8]; // {near, far} x {ll, lr, ur, ul}
	vector<vert_wrap_t> verts;
	get_frustum_corners(pts);

	for (unsigned d = 0; d < 2; ++d) {
		for (unsigned i = 0; i < 4; ++i) {verts.push_back(pts[4*d+i]);} // near and far clipping planes
	}
	for (unsigned i = 0; i < 4; ++i) { // sides
		verts.push_back(pts[0+((i+0)&3)]);
		verts.push_back(pts[0+((i+1)&3)]);
		verts.push_back(pts[4+((i+1)&3)]);
		verts.push_back(pts[4+((i+0)&3)]);
	}
	draw_quad_verts_as_tris(verts);
}


void draw_simple_cube(cube_t const &c, bool texture, unsigned dim_mask, vector3d const *const view_dir) {
	draw_cube(c.get_cube_center(), c.dx(), c.dy(), c.dz(), texture, 1.0, 0, view_dir, dim_mask);
}

// need to do something with tex coords for scale
// Note: cube extends from pos +/- 0.5*(sx, sy, sz)
void draw_cube(point const &pos, float sx, float sy, float sz, bool texture, float texture_scale,
	bool proportional_texture, vector3d const *const view_dir, unsigned dim_mask, bool swap_y_st)
{
	point const scale(sx, sy, sz);
	vector3d const xlate(pos - 0.5*scale); // move origin from center to min corner
	vert_norm_tc verts[24]; // max number of verts that can be drawn is 24, assuming all faces are drawn
	unsigned vix(0);
		
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);
		if (!(dim_mask & (1<<n))) continue;

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (view_dir && (((*view_dir)[n] < 0.0) ^ j)) continue; // back facing
			vector3d norm(zero_vector);
			point pt;
			norm[n] = (2.0*j - 1.0); // -1 or 1
			pt[n]   = j;

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				pt[d[1]] = s1;

				for (unsigned k = 0; k < 2; ++k, ++vix) { // iterate over vertices
					pt[d[0]] = k^j^s1^1; // need to orient the vertices differently for each side
					verts[vix].v = pt*scale + xlate;
					verts[vix].n = norm;

					if (texture) {
						bool const st(swap_y_st ? (i&1) : 0);
						verts[vix].t[ st] = (proportional_texture ? scale[d[1]] : 1.0)*texture_scale*pt[d[1]];
						verts[vix].t[!st] = (proportional_texture ? scale[d[0]] : 1.0)*texture_scale*pt[d[0]];
					}
				} // for k
			} // for s1
		} // for j
	} // for i
	assert(vix <= 24);
	if (vix > 0) {draw_quad_verts_as_tris(verts, vix);}
}

void draw_cube_verts_only(cube_t const &c) { // simplified version of draw_cube() with no normals, texture coordinates, or culling
	point const scale(c.get_size());
	vector3d const xlate(c.get_cube_center() - 0.5*scale); // move origin from center to min corner
	vert_wrap_t verts[24];

	for (unsigned i = 0, vix = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			for (unsigned s1 = 0; s1 < 2; ++s1) {
				point pt;
				pt[n]    = j;
				pt[d[1]] = s1;

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					pt[d[0]] = k^j^s1^1; // need to orient the vertices differently for each side
					verts[vix++] = pt*scale + xlate;
				}
			} // for s1
		} // for j
	} // for i
	draw_quad_verts_as_tris(verts, 24);
}


void gen_quad_tex_coords(float *tdata, unsigned num, unsigned stride) { // stride is in floats

	for (unsigned i = 0, off = 0; i < num; i++) { // 00 10 11 01 for every quad
		tdata[off] = 0.0; tdata[off+1] = 0.0; off += stride;
		tdata[off] = 1.0; tdata[off+1] = 0.0; off += stride;
		tdata[off] = 1.0; tdata[off+1] = 1.0; off += stride;
		tdata[off] = 0.0; tdata[off+1] = 1.0; off += stride;
	}
}


void gen_quad_tri_tex_coords(float *tdata, unsigned num, unsigned stride) { // stride is in floats

	for (unsigned i = 0, off = 0; i < num; i++) { // for every tri
		tdata[off] = -0.5; tdata[off+1] =  1.0; off += stride;
		tdata[off] =  0.5; tdata[off+1] = -1.0; off += stride;
		tdata[off] =  1.5; tdata[off+1] =  1.0; off += stride;
	}
}


// ******************** Sphere VBOs ********************


void free_sphere_vbos() {
	delete_and_zero_vbo(predef_sphere_vbo);
}


void setup_sphere_vbos() {

	assert(MAX_SPHERE_VBO_NDIV > 0);
	if (predef_sphere_vbo > 0) return; // already finished
	vector<sd_sphere_d::vertex_type_t> verts;
	sphere_point_norm spn;
	sphere_vbo_offsets[0] = 0;

	for (unsigned i = 3; i <= MAX_SPHERE_VBO_NDIV; ++i) {
		for (unsigned half = 0; half < 2; ++half) {
			for (unsigned tex = 0; tex < 2; ++tex) {
				sd_sphere_d sd(all_zeros, 1.0, i);
				sd.gen_points_norms(spn, 0.0, 1.0, 0.0, (half ? 0.5 : 1.0));
				sd.get_triangle_strip_pow2(verts, 1, 0.0, 1.0, 0.0, (half ? 0.5 : 1.0));
				sphere_vbo_offsets[((i-1) << 2) + (half << 1) + tex] = verts.size();
			}
		}
	}
	create_vbo_and_upload(predef_sphere_vbo, verts, 0, 1); // ~8MB
}


void bind_draw_sphere_vbo(bool textured, bool normals) {

	assert(predef_sphere_vbo > 0);
	bind_vbo(predef_sphere_vbo);
	set_array_client_state(1, textured, normals, 0);
	sd_sphere_d::vertex_type_t::set_vbo_arrays(0);
}


void draw_sphere_vbo_pre_bound(int ndiv, bool textured, bool half, unsigned num_instances) {

	assert(ndiv >= 3);
	assert(ndiv <= (int)MAX_SPHERE_VBO_NDIV);
	unsigned const ix(((ndiv-1) << 2) + (half << 1) + textured), off1(sphere_vbo_offsets[ix-1]), off2(sphere_vbo_offsets[ix]);
	assert(off1 < off2);
	check_mvm_update();
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, off1, (off2 - off1), num_instances); // uses triangle strips separated by degenerate triangles
	++num_frame_draw_calls;
}


bool sphere_vbo_bound(0);

void begin_sphere_draw(bool textured) {
	assert(!sphere_vbo_bound); sphere_vbo_bound = 1;
	bind_draw_sphere_vbo(textured, 1);
}
void end_sphere_draw() {
	assert(sphere_vbo_bound); sphere_vbo_bound = 0;
	bind_vbo(0);
}
void draw_sphere_vbo_raw(int ndiv, bool textured, bool half, unsigned num_instances) {
	if (!sphere_vbo_bound) {bind_draw_sphere_vbo(textured, 1);}
	draw_sphere_vbo_pre_bound(ndiv, textured, half, num_instances);
	if (!sphere_vbo_bound) {bind_vbo(0);}
}


void draw_sphere_vbo(point const &pos, float radius, int ndiv, bool textured, bool half, bool bfc) {

	if (ndiv <= (int)MAX_SPHERE_VBO_NDIV) { // speedup is highly variable
		assert(ndiv > 0);
		bool const has_xform(radius != 1.0 || pos != all_zeros);
		//if (shader_loc >= 0) {shader_t::set_uniform_vector4d(shader_loc, vector4d(pos, radius));}
		
		if (has_xform) {
			fgPushMatrix();
			translate_to(pos);
			uniform_scale(radius);
		}
		draw_sphere_vbo_raw(ndiv, textured, half);
		if (has_xform) {fgPopMatrix();}
	}
	else if (half) {draw_subdiv_sphere_section(pos, radius, ndiv, textured, 0.0, 1.0, 0.0, 0.5);}
	else {draw_subdiv_sphere(pos, radius, ndiv, textured, !bfc);}
}


void draw_sphere_vbo_back_to_front(point const &pos, float radius, int ndiv, bool textured, bool enable_front, bool enable_back) {

	glEnable(GL_CULL_FACE);

	if (ndiv <= (int)MAX_SPHERE_VBO_NDIV && (radius != 1.0 || pos != all_zeros)) { // optimization for common case (shared transforms and VBO bind)
		bind_draw_sphere_vbo(textured, 1);
		fgPushMatrix();
		translate_to(pos);
		uniform_scale(radius);

		for (unsigned i = unsigned(!enable_back); i < unsigned(1+enable_front); ++i) { // kind of slow
			glCullFace(i ? GL_BACK : GL_FRONT);
			draw_sphere_vbo_pre_bound(ndiv, textured);
		}
		fgPopMatrix();
		bind_vbo(0);
	}
	else {
		for (unsigned i = unsigned(!enable_back); i < unsigned(1+enable_front); ++i) { // kind of slow
			glCullFace(i ? GL_BACK : GL_FRONT);
			draw_sphere_vbo(pos, radius, ndiv, textured); // cull?, partial sphere?
		}
	}
	glDisable(GL_CULL_FACE);
}

