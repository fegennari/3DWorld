// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "mesh.h"
#include "transform_obj.h"
#include "physics_objects.h"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>


extern float base_gravity, tstep, fticks;
extern obj_type object_types[];


// *** xform_matrix ***


float       *xform_matrix::get_ptr()       {glm::mat4       &m(*this); return glm::value_ptr(m);}
float const *xform_matrix::get_ptr() const {glm::mat4 const &m(*this); return glm::value_ptr(m);}


void xform_matrix::normalize() {

	float *m(get_ptr());

	for (unsigned i = 0; i < 3; ++i) { // renormalize matrix to account for fp error
		float const dist(sqrt(m[i+0]*m[i+0] + m[i+4]*m[i+4] + m[i+8]*m[i+8]));
		m[i+0] /= dist;
		m[i+4] /= dist;
		m[i+8] /= dist;
	}
}


matrix_stack_t mvm_stack; // modelview matrix
matrix_stack_t pjm_stack; // projection matrix
int matrix_mode(0); // 0 = MVM, 1 = PJM
bool mvm_changed(0);

matrix_stack_t &get_matrix_stack() {return (matrix_mode ? pjm_stack : mvm_stack);}

void mark_matrix_changed() {
	if (!matrix_mode) {mvm_changed = 1;}
}

void fgMatrixMode(int val) {
	if      (val == GL_PROJECTION) {matrix_mode = 1;}
	else if (val == GL_MODELVIEW ) {matrix_mode = 0;}
	else {assert(0);}
}

void fgPushMatrix  () {get_matrix_stack().push();} // matrix not change
void fgPopMatrix   () {get_matrix_stack().pop();      mark_matrix_changed();}
void fgLoadIdentity() {get_matrix_stack().identity(); mark_matrix_changed();}

void fgTranslate(float x, float y, float z) {
	get_matrix_stack().assign(glm::translate(get_matrix_stack().top(), glm::vec3(x, y, z)));
	mark_matrix_changed();
}
void fgScale(float x, float y, float z) {
	get_matrix_stack().assign(glm::scale(get_matrix_stack().top(), glm::vec3(x, y, z)));
	mark_matrix_changed();
}
void fgScale(float s) {fgScale(s, s, s);}

void fgRotate(float angle, float x, float y, float z) {
	get_matrix_stack().assign(glm::rotate(get_matrix_stack().top(), angle, glm::vec3(x, y, z)));
	mark_matrix_changed();
}
void fgRotateDegrees(float angle, float x, float y, float z) {fgRotate(TO_RADIANS*angle, x, y, z);}

void fgPerspective(float fov_y, float aspect, float near_clip, float far_clip) {
	get_matrix_stack().assign(glm::perspective(TO_RADIANS*fov_y, aspect, near_clip, far_clip));
	mark_matrix_changed();
}

void fgLookAt(float eyex, float eyey, float eyez, float centerx, float centery, float centerz, float upx, float upy, float upz) {
	get_matrix_stack().assign(glm::lookAt(glm::vec3(eyex, eyey, eyez), glm::vec3(centerx, centery, centerz), glm::vec3(upx, upy, upz)));
	mark_matrix_changed();
}

void fgMultMatrix(xform_matrix const &m) {
	get_matrix_stack().assign(m * get_matrix_stack().top());
}

xform_matrix fgGetMVM() {return mvm_stack.top();}
xform_matrix fgGetPJM() {return pjm_stack.top();}

#ifdef USE_FG_TRANSFORMS
void xform_matrix::assign_mv_from_gl() {*this = fgGetMVM();}
void xform_matrix::assign_pj_from_gl() {*this = fgGetPJM();}
#else
void xform_matrix::assign_mv_from_gl() {glGetFloatv(GL_MODELVIEW_MATRIX,  get_ptr());}
void xform_matrix::assign_pj_from_gl() {glGetFloatv(GL_PROJECTION_MATRIX, get_ptr());}
#endif


// *** mesh2d ***


void mesh2d::clear() {
	
	delete [] pmap;
	delete [] rmap;
	delete [] ptsh;
	pmap = NULL;
	rmap = NULL;
	ptsh = NULL;
	size = 0;
}


void mesh2d::set_size(unsigned sz) {

	assert(sz > 0);
	assert(size == 0 || sz == size);
	clear();
	size = sz;
}


template<typename T> void mesh2d::alloc_ptr(T *&p, T const val) {

	unsigned const num(get_num());
	delete [] p;
	p = new T[num];
	for (unsigned i = 0; i < num; ++i) p[i] = val;
}


void mesh2d::alloc_pmap() {alloc_ptr(pmap, 0.0f);}
void mesh2d::alloc_rmap() {alloc_ptr(rmap, bool(1));}
void mesh2d::alloc_emap() {alloc_ptr(emap, 0.0f);}
void mesh2d::alloc_ptsh() {alloc_ptr(ptsh, all_zeros);}


void mesh2d::reset_pmap() {

	if (!pmap) {alloc_pmap(); return;} // will be reset
	unsigned const num(get_num());
	for (unsigned i = 0; i < num; ++i) pmap[i] = 0.0;
}


void mesh2d::add_random(float mag, float min_mag, float max_mag, unsigned skipval) {

	if (!pmap) alloc_pmap();
	unsigned const num(get_num());

	for (unsigned i = (rand()%(skipval+1)); i < num; i += (skipval+1)) {
		pmap[i] = max(min_mag, min(max_mag, (pmap[i] + mag*signed_rand_float())));
	}
}


void mesh2d::mult_by(float val) {

	if (!pmap) return;
	unsigned const num(get_num());
	for (unsigned i = 0; i < num; ++i) pmap[i] *= val;
}


void mesh2d::unset_rand_rmap(unsigned num_remove) {

	if (!rmap) alloc_rmap();
	for (unsigned i = 0; i < num_remove; ++i) rmap[choose_rand()] = 0; // doesn't check for already removed elements
}


void mesh2d::set_rand_expand(float mag, unsigned num_exp) {

	if (!emap) alloc_emap();
	for (unsigned i = 0; i < num_exp; ++i) emap[choose_rand()] += mag; // doesn't check for already removed elements
}


void mesh2d::set_rand_translate(point const &tp, unsigned num_trans) {

	if (tp == all_zeros) return;
	if (!ptsh) alloc_ptsh();
	for (unsigned i = 0; i < num_trans; ++i) ptsh[choose_rand()] += tp; // doesn't check for already translated elements
}


void mesh2d::draw_perturbed_sphere(point const &pos, float radius, int ndiv, bool tex_coord) const {

	if (!pmap && !rmap && !emap && !ptsh && expand == 0.0) {
		draw_sphere_vbo(pos, radius, ndiv, 1);
	}
	else { // ndiv unused
		if (pmap || rmap || ptsh || emap) assert(size > 0);
		point const camera(get_camera_all());
		draw_subdiv_sphere(pos, radius, size, camera, pmap, tex_coord, 1, rmap, emap, ptsh, expand);
	}
}


// *** transform_data ***


void transform_data::set_perturb_size(unsigned i, unsigned sz) {

	assert(i < perturb_maps.size());

	if (perturb_maps[i].pmap) {
		assert(perturb_maps[i].get_size() == sz);
	}
	else {
		perturb_maps[i].set_size(sz);
		perturb_maps[i].alloc_pmap();
	}
}


void transform_data::add_perturb_at(unsigned s, unsigned t, unsigned i, float val, float min_mag, float max_mag) {

	assert(i < perturb_maps.size());
	perturb_maps[i].set_val(s, t, max(min_mag, min(max_mag, (perturb_maps[i].get_val(s, t) + val))));
}


void transform_data::reset_perturb_if_set(unsigned i) {
	
	assert(i < perturb_maps.size());
	if (perturb_maps[i].pmap) perturb_maps[i].reset_pmap();
}


transform_data::~transform_data() {
		
	for (unsigned i = 0; i < perturb_maps.size(); ++i) {perturb_maps[i].clear();}
}


// *** deformation code ***


void apply_obj_mesh_roll(xform_matrix &matrix, shader_t &shader, point const &pos, point const &lpos, float radius, float a_add, float a_mult) {

	if (pos != lpos) {
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

		if (!point_outside_mesh(xpos, ypos)) {
			vector3d const delta(pos, lpos);
			float const dmag(delta.mag()), angle(a_mult*(dmag/radius) + a_add);
			vector3d const vrot(cross_product(surface_normals[ypos][xpos], delta/dmag));

			if (vrot.mag() > TOLERANCE) {
				matrix.normalize();
				matrix = glm::rotate(glm::mat4(), angle, vec3_from_vector3d(vrot)) * matrix;
			}
		}
	}
	// TODO: set some shader uniform from this matrix
#ifdef USE_FG_TRANSFORMS
	fgMultMatrix(matrix);
#else
	glMultMatrixf(matrix.get_ptr());
#endif
}


void deform_obj(dwobject &obj, vector3d const &norm, vector3d const &v0) { // apply collision deformations

	float const deform(object_types[obj.type].deform);
	if (deform == 0.0) return;
	assert(deform > 0.0 && deform < 1.0);
	vector3d const vd(obj.velocity, v0);
	float const vthresh(base_gravity*GRAVITY*tstep*object_types[obj.type].gravity), vd_mag(vd.mag());

	if (vd_mag > max(2.0f*vthresh, 12.0f/fticks) && (fabs(v0.x) + fabs(v0.y)) > 0.01) { // what about when it hits the ground/mesh?
		float const deform_mag(SQRT3*deform*min(1.0, 0.05*vd_mag));
		UNROLL_3X(obj.vdeform[i_] -= fabs(norm[i_])*deform_mag;)
		obj.vdeform *= SQRT3/obj.vdeform.mag(); // normalize the volume
		float vdmin(1.0 - deform);
		UNROLL_3X(vdmin = min(vdmin, obj.vdeform[i_]);)

		if (vdmin < (1.0 - deform)) {
			UNROLL_3X(obj.vdeform[i_] += ((1.0 - deform) - vdmin);)
			obj.vdeform *= SQRT3/obj.vdeform.mag(); // re-normalize
		}
	}
}


void update_deformation(dwobject &obj) {

	if (obj.vdeform != all_ones && object_types[obj.type].def_recover > 0.0) {
		obj.vdeform += all_ones*(fticks*object_types[obj.type].def_recover);
		obj.vdeform *= SQRT3/obj.vdeform.mag(); // normalize the volume
	}
}


void mirror_about_plane(vector3d const &norm, point const &pt) { // applies to GL state

	float const dp(dot_product(pt, norm));
	float const m[16] = {1-2*norm.x*norm.x,  -2*norm.x*norm.y,  -2*norm.x*norm.z, 0.0,
			              -2*norm.x*norm.y, 1-2*norm.y*norm.y,  -2*norm.y*norm.z, 0.0,
		                  -2*norm.x*norm.z,  -2*norm.y*norm.z, 1-2*norm.z*norm.z, 0.0,
		                   2*dp*norm.x,       2*dp*norm.y,       2*dp*norm.z,     1.0};
#ifdef USE_FG_TRANSFORMS
	fgMultMatrix(glm::make_mat4(m));
#else
	glMultMatrixf(m);
#endif
}


// **************** INSTANCING ****************


void instance_render_t::add_cur_inst() {

	xform_matrix xf;
	xf.assign_mv_from_gl();
	add_inst(xf);
}

// Note: instance_render_t::draw_and_clear() is defined in shaders.cpp


