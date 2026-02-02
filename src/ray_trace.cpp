// 3D World - Ray Tracing Code
// by Frank Gennari
// 2/14/10
#include "3DWorld.h"
#include "lightmap.h"
#include "mesh.h"
#include "model3d.h"
#include "binary_file_io.h"
#include <atomic>
#include <thread>
#include <omp.h>


bool const COLOR_FROM_COBJ_TEX = 0; // 0 = fast/average color, 1 = true color
float const RAY_WEIGHT    = 4.0E5;
float const WEIGHT_THRESH = 0.01;
float const DIFFUSE_REFL  = 0.9; // 90%  diffuse  reflectivity
float const SPEC_REFL     = 1.0; // 100% specular reflectivity
float const SNOW_ALBEDO   = 0.9;
float const ICE_ALBEDO    = 0.8;

bool keep_beams(0); // debugging mode
bool kill_raytrace(0);
bool no_stat_moving(0); // generally not thread safe for dynamic lighting update, since BVH is rebuilt per-frame; also, wrong to cache lighting for moving cobjs
unsigned NPTS(50000), NRAYS(40000), LOCAL_RAYS(1000000), GLOBAL_RAYS(1000000), DYNAMIC_RAYS(1000000), NUM_THREADS(1), MAX_RAY_BOUNCES(20);
std::atomic<unsigned long long> tot_rays(0), num_hits(0), cells_touched(0);
unsigned const NUM_RAY_SPLITS [NUM_LIGHTING_TYPES] = {1, 1, 1, 1, 1}; // sky, global, local, cobj_accum, dynamic
unsigned const INIT_RAY_SPLITS[NUM_LIGHTING_TYPES] = {1, 4, 1, 1, 1}; // sky, global, local, cobj_accum, dynamic

extern bool has_snow, combined_gu, global_lighting_update, lighting_update_offline, store_cobj_accum_lighting_as_blocked;
extern int read_light_files[], write_light_files[], display_mode, DISABLE_WATER;
extern float water_plane_z, temperature, snow_depth, ray_step_size_mult, first_ray_weight[];
extern char *lighting_file[];
extern point sun_pos, moon_pos;
extern vector<light_source> light_sources_a;
extern vector<light_source_trig> light_sources_d;
extern coll_obj_group coll_objects;
extern platform_cont platforms;
extern vector<beam3d> beams;
extern lmap_manager_t lmap_manager;
extern llv_vect local_light_volumes;
extern indir_dlight_group_manager_t indir_dlight_group_manager;
extern cube_light_src_vect sky_cube_lights, global_cube_lights;
extern model3ds all_models;


struct face_ray_accum_t {

	colorRGBA color; // +weight
	cube_t bcube;
	unsigned num_rays;

	face_ray_accum_t() : color(0,0,0,0), bcube(all_zeros), num_rays(0) {}

	void add_ray(point const &pos, colorRGBA const &c, float weight) {
		if (num_rays == 0) {bcube.set_from_point(pos);} else {bcube.union_with_pt(pos);}
		color += c * weight;
		++num_rays;
	}
	void add(face_ray_accum_t const &v) {
		if (v.num_rays == 0) return; // nothing to add
		if (num_rays == 0) {bcube = v.bcube;} else {bcube.union_with_cube(v.bcube);}
		color += v.color;
		num_rays += v.num_rays;
	}
};

struct rt_ray_t { // size = 24 (40 if uncompressed)
	point pos;
	norm_comp dir;
	color_wrapper color;
	float weight;

	rt_ray_t() : weight(0.0) {}
	rt_ray_t(point const &p1, point const &p2, colorRGB color_, float weight_) : pos(p1), weight(weight_) {
		dir.set_norm((p2 - p1).get_norm());
		float const max_comp(color_.get_max_component());
		color_ *= 1.0/max_comp;
		weight *= max_comp;
		color.set_c4(color_);
	}
	point get_p2(float length) const {return pos + length*dir.get_norm();}
	colorRGBA get_color() const {return color.get_c4();}
};

struct cobj_ray_accum_t {

	face_ray_accum_t vals[6]; // one per cube face
	vector<rt_ray_t> rays;

	void add_ray(point const &p1, point const &p2, colorRGB const &color, float weight, unsigned face) {
		assert(face < 6);
		vals[face].add_ray(p1, color, weight);
		rays.push_back(rt_ray_t(p1, p2, color, weight));
	}
	unsigned get_count() const {
		unsigned n(0);
		for (unsigned i = 0; i < 6; ++i) {n += vals[i].num_rays;}
		return n;
	}
	cube_t get_bcube(bool const expand=0) const {
		cube_t bcube(all_zeros);
		bool is_first(0);
		for (unsigned i = 0; i < 6; ++i) {
			if (vals[i].num_rays == 0) continue;
			if (is_first) {bcube = vals[i].bcube; is_first = 0;} else {bcube.union_with_cube(vals[i].bcube);}
		}
		if (expand) {bcube.expand_by(1.0E-6);} // expand slightly to turn adjacency into intersection
		return bcube;
	}
	void merge(cobj_ray_accum_t const &v) {
		for (unsigned i = 0; i < 6; ++i) {vals[i].add(v.vals[i]);}
		rays.insert(rays.end(), v.rays.begin(), v.rays.end());
	}
};

unsigned const magic_val = 0xbeefdead;

struct cobj_ray_accum_map_t : public map<unsigned, cobj_ray_accum_t> {

	void add_ray(unsigned id, point const &p1, point const &p2, colorRGB const &color, float weight, unsigned face) {
		operator[](id).add_ray(p1, p2, color, weight, face);
	}
	void merge(cobj_ray_accum_map_t const &m) { // merge maps across threads
		for (const_iterator i = m.begin(); i != m.end(); ++i) {operator[](i->first).merge(i->second);}
	}
	bool read(FILE *fp) {
		clear();
		unsigned sz(0), magic(0);
		if (fread(&magic, sizeof(unsigned), 1, fp) != 1) return 0; // read number of entries
		if (magic != magic_val) {cerr << "Incorrect cobj ray accumulation file type" << endl; return 0;}
		if (fread(&sz, sizeof(unsigned), 1, fp) != 1) return 0; // read number of entries
		for (unsigned i = 0; i < sz; ++i) {
			pair<unsigned, cobj_ray_accum_t> val;
			if (fread(&val.first,  sizeof(unsigned),         1, fp) != 1) return 0; // read ID
			if (fread(&val.second, sizeof(face_ray_accum_t), 6, fp) != 6) return 0; // read 6 values
			unsigned nrays(0);
			if (fread(&nrays, sizeof(unsigned), 1, fp) != 1) return 0; // read number of rays
			val.second.rays.resize(nrays);
			if (fread(&val.second.rays.front(), sizeof(rt_ray_t), nrays, fp) != nrays) return 0; // read ray data
			bool const did_ins(insert(val).second);
			assert(did_ins); // no duplicate ids
		}
		return 1;
	}
	bool write(FILE *fp) const {
		if (fwrite(&magic_val, sizeof(unsigned), 1, fp) != 1) return 0; // write magic value
		unsigned const sz(size());
		if (fwrite(&sz, sizeof(unsigned), 1, fp) != 1) return 0; // write number of entries
		for (const_iterator i = begin(); i != end(); ++i) {
			if (fwrite(&i->first,  sizeof(unsigned),         1, fp) != 1) return 0; // write ID
			if (fwrite(&i->second, sizeof(face_ray_accum_t), 6, fp) != 6) return 0; // write 6 values
			unsigned const nrays(i->second.rays.size());
			if (fwrite(&nrays, sizeof(unsigned), 1, fp) != 1) return 0; // write number of rays
			if (fwrite(&i->second.rays.front(), sizeof(rt_ray_t), nrays, fp) != nrays) return 0; // write ray data
		}
		return 1;
	}
	void open_and_read(string const &filename, bool show_stats=0) {
		FILE *fp(fopen(filename.c_str(), "rb")); // read a binary file
		if (fp == nullptr) {cerr << "failed to open file '" << filename << "' for reading" << endl; exit(1);}
		bool const success(read(fp));
		checked_fclose(fp);
		if (!success) {cerr << "failed to read file '" << filename << "'" << endl; exit(1);}
		cout << "Read cobj accum lighting file " << filename << endl;
		if (show_stats) {stats();}
	}
	void open_and_write(string const &filename, bool show_stats=0) const {
		FILE *fp(fopen(filename.c_str(), "wb")); // write a binary file
		if (fp == nullptr) {cerr << "failed to open file '" << filename << "' for writing" << endl; exit(1);}
		bool const success(write(fp));
		checked_fclose(fp);
		if (!success) {cerr << "failed to write file '" << filename << "'" << endl; exit(1);}
		cout << "Wrote cobj accum lighting file " << filename << endl;
		if (show_stats) {stats();}
	}
	void stats() const {
		cout << "cobj lighting stats:" << endl;
		for (auto i = begin(); i != end(); ++i) {
			cout << "cobj: " << i->first << ", count: " << i->second.get_count() << ", rays stored: " << i->second.rays.size() << endl;
			for (unsigned n = 0; n < 6; ++n) {
				face_ray_accum_t const &val(i->second.vals[n]);
				if (val.num_rays == 0) continue; // no rays for this face
				cout << "dim: " << (n>>1) << ", dir: " << (n&1) << ", num_rays: " << val.num_rays << endl;
				cout << "color: " << val.color.str() << endl;
				cout << "bcube: " << val.bcube.str() << endl;
			}
		}
	}
};

cobj_ray_accum_map_t merged_accum_map;


float get_scene_radius() {return sqrt(2.0f*(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE + Z_SCENE_SIZE*Z_SCENE_SIZE));}
float get_step_size()    {return 0.3f*ray_step_size_mult*(DX_VAL + DY_VAL + DZ_VAL);}

void increment_printed_number(unsigned num) {

	for (unsigned n = max(num, 1U); n > 0; n /= 10) {cout << "\b";}
	cout << (num+1);
	cout.flush();
}

bool is_ltype_dynamic(int ltype) {return (ltype >= LIGHTING_DYNAMIC);}
int clamp_ltype_range(int ltype) {return (is_ltype_dynamic(ltype) ? LIGHTING_DYNAMIC : ltype);}
bool enable_platform_lights(int ltype) {return (ltype == LIGHTING_SKY);} // sky lighting only for now

light_volume_local &get_local_light_volume(int ltype) {

	assert(is_ltype_dynamic(ltype)); // it's a local lighting volume
	unsigned const llvol_ix(ltype - LIGHTING_DYNAMIC);
	assert(llvol_ix < local_light_volumes.size());
	return *local_light_volumes[llvol_ix];
}


// Note: weight can be negative
unsigned add_path_to_lmcs(lmap_manager_t *lmgr, cube_t *bcube, point p1, point const &p2, float weight, colorRGBA const &color, int ltype, bool first_pt) {

	bool const dynamic(is_ltype_dynamic(ltype));
	if (first_pt && dynamic) return 0; // since dynamic lights already have a direct lighting component, we skip the first ray here to avoid double counting it
	if (first_pt) {weight *= first_ray_weight[ltype];} // lower weight - handled by direct illumination
	if (fabs(weight) < TOLERANCE) return 0;
	weight *= ray_step_size_mult;
	colorRGBA const cw(color*weight);
	unsigned const nsteps(1 + unsigned(p2p_dist(p1, p2)/get_step_size())); // round up (dist can be 0)
	vector3d const step((p2 - p1)/nsteps); // at least two points
	if (!first_pt) {p1 += step;} // move past the first step so we don't double count

	// better time vs. quality tradeoff using a proper line drawing algorithm that chooses step size by distance to closest grid boundary and multiplies weight by segment length?
	if (dynamic) { // it's a local lighting volume
		light_volume_local &lvol(get_local_light_volume(ltype));

		for (unsigned s = 0; s < nsteps; ++s) {
			lvol.add_color(p1, cw);
			p1 += step;
		}
	}
	else { // use the lmgr
		assert(lmgr != nullptr && lmgr->is_allocated());

		for (unsigned s = 0; s < nsteps; ++s) {
			lmcell *lmc(lmgr->get_lmcell_round_down(p1));
		
			if (lmc != NULL) { // could use a mutex here, but it seems too slow
				float *color(lmc->get_offset(ltype));
				ADD_LIGHT_CONTRIB(cw, color);
				if (ltype != LIGHTING_LOCAL) {color[3] += weight;}
			}
			p1 += step;
		}
		if (bcube) {
			bcube->assign_or_union_with_pt(p1);
			bcube->union_with_pt(p2);
		}
		lmgr->was_updated = 1;
	}
	return nsteps;
}


void cast_light_ray(lmap_manager_t *lmgr, point p1, point p2, float weight, float weight0, colorRGBA color, float line_length,
	int ignore_cobj, int ltype, unsigned depth, rand_gen_t &rgen, cobj_ray_accum_map_t *accum_map, cube_t *bcube=nullptr)
{
	if (depth > MAX_RAY_BOUNCES) return;
	if (ltype == LIGHTING_DYNAMIC && depth > 4) return; // use a sensible default since this is running during rendering
	//assert(!is_nan(p1) && !is_nan(p2));
	++tot_rays;

	// find intersection point with scene cobjs
	point orig_p1(p1);
	if (!do_line_clip_scene(p1, p2, min(zbottom, czmin), max(ztop, czmax))) return;
	if ((display_mode & 0x01) && is_under_mesh(p1)) return;
	int cindex(-1), xpos(0), ypos(0);
	point cpos(p2);
	vector3d cnorm;
	float t(0.0), zval(0.0);
	bool snow_coll(0), ice_coll(0), water_coll(0), mesh_coll(0);
	vector3d const dir((p2 - p1).get_norm());
	bool coll(check_coll_line_exact(p1, p2, cpos, cnorm, cindex, 0.0, ignore_cobj, 1, 0, 1, 1, (p1 == orig_p1), no_stat_moving)); // fast=1, exclude voxels, maybe skip init colls
	assert(coll ? (cindex >= 0 && cindex < (int)coll_objects.size()) : (cindex == -1));

	// find the intersection point with the model3ds
	colorRGBA model_color;
	bool const model_coll(all_models.check_coll_line(p1, cpos, cpos, cnorm, model_color, 1));
	coll |= model_coll;

	// find intersection point with mesh (approximate)
	// Note: the !coll test is a big optimization but not entirely correct, as we can have a ray that intersects the mesh and then hits a cobj below the mesh;
	// however, if the scene is properly built, the cobj surface should not be visible anyway
	if ((display_mode & 0x01) && !coll && p1.z != p2.z && line_intersect_mesh(p1, p2, xpos, ypos, zval, 0, 0)) {
		assert(!point_outside_mesh(xpos, ypos));
		
		if (!is_mesh_disabled(xpos, ypos)) {
			if (p2.z >= p1.z) return; // starts under mesh = bad
			cpos  = (p1 + (p2 - p1)*(zval + SMALL_NUMBER - p1.z)/(p2.z - p1.z));
			cnorm = vertex_normals[ypos][xpos];
			coll  = mesh_coll = 1;
		}
	}

	// intersection with water/ice
	if (!DISABLE_WATER && coll && p1.z >= water_plane_z && cpos.z < water_plane_z) { // what if p1.z < water_plane_z?
		if (temperature <= W_FREEZE_POINT) { // ice
			float const t((water_plane_z - cpos.z)/(p1.z - cpos.z));
			cpos    += (p1 - cpos)*t;
			cnorm    = plus_z; // always up? what about ripples in ice?
			ice_coll = 1;
		}
		else { // water
			water_coll = 1;
		}
	}

	// intersection with snow (not exactly correct for curved surfaces)
	if (coll && has_snow && cnorm.z > 0.0) {
		float zval;
		vector3d snow_cnorm(cnorm);

		// could iterate, but cpos won't necessarily converge because snow strips are discrete rather than smooth/interpolated
		if (get_snow_height(cpos, 2.0*snow_depth, zval, snow_cnorm, 0) && zval >= cpos.z) {
			assert(snow_cnorm.z >= 0.0);
			vector3d const delta(p1 - cpos), delta_dir(delta.get_norm());
			float const height(zval - cpos.z), dp(dot_product(snow_cnorm, delta_dir));
			
			if (dp > 0.0) { // else ray entered under the snow, so we can't really update the position???
				float const dist(min(delta.mag(), height*snow_cnorm.z/dp));
				snow_coll = 1;
				cpos     += delta_dir*dist;
				cnorm     = snow_cnorm; // should we do this?
			}
		}
	}
	point p_end(p2);
	if ( coll) {p2 = cpos;}
	if (keep_beams && p1 != p2) {beams.push_back(beam3d(!coll, 1, p1, p2, color, 0.1*weight));} // TESTING
	if (!coll) return; // more efficient to do this up here and let a reverse ray from the sky light this path

	// walk from p1 to p2, adding light to all lightmap cells encountered
	cells_touched += add_path_to_lmcs(lmgr, bcube, p1, p2, weight, color, ltype, (depth == 0));
	++num_hits;
	//if (!coll)    return;
	if (p1 == p2) return; // line must have started inside a cobj - this is bad, but what can we do?

	// update ray weight based on object material properties
	float specular(0.0), shine(1.0);

	if (water_coll) {
		assert(!ice_coll && !snow_coll);
		// We want full Fresnel reflected/refracted rays with optical path attenuation through water
		// such as in water_surface_draw::blend_reflection_color()
		// However, we can make simplifying assumptions that make the code much simpler and more efficient:
		// 1. Most of the light is refracted though clear water, so we can ignore the ~5% reflected light
		// 2. The water is generally pretty shallow in most scenes and will reflect diffusely on the bottom,
		//    so we can ignore the angle of refraction (and therefore skip the intersection point entirely)
		// 3. This is a preproc step, so we can ignore dynamic water and impurities that haven't been added yet
		// 4. We can assume the optical path of the reflected ray is similar to the optical path of the incident ray
		// 5. The water is clear, so we can ignore scattering effects
		// Therefore, we simply multiply by the base water color and attenuate by 2x the incident optical path
		vector3d const delta(p2 - p1);
		if (delta.z > -TOLERANCE) return; // too shallow of an angle, assume attenuated to nothing
		float const wdepth(water_plane_z - cpos.z);
		float const dist(-2*wdepth*delta.mag()/delta.z); // multiply by 2 to account for both directions
		assert(dist >= 0.0);
		colorRGBA water_color(WATER_C);
		water_color.A = 1.0;  // make solid for volume attenuation (not surface)
		water_color  *= 0.95; // refracted component
		atten_by_water_depth(&(water_color.R), 0.8*dist); // assumes dist scaling for outside water
		weight *= water_color.get_luminance();
		color   = color.modulate_with(water_color);
	} // Note: no else
	if (snow_coll) {
		weight  *= SNOW_ALBEDO; // white
		specular = 0.5;
		shine    = 50.0;
	}
	else if (ice_coll) {
		weight  *= ICE_ALBEDO*ICE_C.get_luminance();
		color    = color.modulate_with(ICE_C);
		specular = 0.5; // see w_spec
		shine    = 60.0;
	}
	else if (mesh_coll) { // collision with mesh
		colorRGBA const lc(get_landscape_texture_color(xpos, ypos));
		weight *= DIFFUSE_REFL*lc.get_luminance(); // 90% diffuse reflectivity
		color   = color.modulate_with(lc);
	}
	else if (model_coll) {
		weight  *= model_color.get_luminance();
		color    = color.modulate_with(model_color); // Note: model texture coords not available here - use average color
		// TODO: get this from the model?
		// requires storing material id in each coll_tquad and calculating approx specular component from specular color and textures
		specular = 0.0;
		shine    = 1.0;
	}
	else { // collision with cobj
		assert(cindex >= 0);
		coll_obj const &cobj(coll_objects[cindex]);

		if (accum_map && enable_platform_lights(ltype) && cobj.is_update_light_platform()) {
			assert(cobj.type == COLL_CUBE); // not yet supported
			assert(!cobj.is_semi_trans()); // not yet supported
			if (fabs(weight) < 2.0*WEIGHT_THRESH*weight0) return; // higher thresh for cobj rays
			unsigned const dim(get_max_dim(cnorm)); // cube intersection dim
			bool const dir(cnorm[dim] > 0); // cube intersection dir
			unsigned const face((dim<<1) + dir);
			accum_map->add_ray(cindex, cpos, p_end, color, weight, face); // use pre-reflection color and weight
			return; // terminate the ray here - will be continued from the hit surface in a later dynamic pass
		}
		else if (cobj.platform_id >= 0 || cobj.is_movable()) {
			// this cobj isn't static, so maybe shouldn't be included - do we skip it in the line intersection check?
			// maybe it's okay to keep this cobj, since it represents the intial object states, and will be correct if it's never moved
		}
		colorRGBA const cobj_color(cobj.get_color_at_point(cpos, cnorm, !COLOR_FROM_COBJ_TEX));
		float const alpha(cobj_color.alpha);
		specular = cobj.cp.spec_color.get_luminance();
		shine    = cobj.cp.shine;
		weight  *= cobj_color.get_luminance();
		color    = color.modulate_with(cobj_color);

		// calculate refracted ray
		if (alpha < 1.0) { // semi-transparent
			float rweight(alpha);

			if (cobj.cp.refract_ix != 1.0) {
				rweight = get_reflected_weight(get_fresnel_reflection(dir, -cnorm, 1.0, cobj.cp.refract_ix), alpha);
			}
			float tweight((1.0 - rweight)*weight); // refracted weight
			
			if (fabs(tweight) > WEIGHT_THRESH*weight0) {
				bool no_transmit(0);

				if (cobj.cp.refract_ix != 1.0) { // refracted
					vector3d v_refract, v_refract2;
					
					if (calc_refraction_angle(dir, v_refract, cnorm, 1.0, cobj.cp.refract_ix)) {
						point const enter_pt(p2);
						p_end = (p2 + v_refract*line_length);
						vector3d cnorm2;

						// test for collision with reversed ray to get the other intersection point
						if (cobj.line_int_exact(p_end, p2, t, cnorm2)) { // not sure what to do if fails or tmax >= 1.0
							point const p_int(p_end + (p2 - p_end)*t);

							if (!dist_less_than(p2, p_int, get_step_size())) {	
								cells_touched += add_path_to_lmcs(lmgr, bcube, p2, p_int, weight, color, ltype, (depth == 0));
								++num_hits;
							}
							if (calc_refraction_angle(v_refract, v_refract2, -cnorm2, cobj.cp.refract_ix, 1.0)) {
								p2    = p_int;
								p_end = p2 + v_refract2*line_length;
								tweight    *= cobj.get_light_transmit(enter_pt, p_int); // can we use p2p_dist(enter_pt, p_int) directly?
								no_transmit = !(fabs(tweight) > WEIGHT_THRESH*weight0);
							}
							else {
								no_transmit = 1; // total internal reflection (could process an internal reflection)
							}
						}
					}
					else {
						no_transmit = 1; // total internal reflection (could process an internal reflection)
					}
				}
				if (!no_transmit) {cast_light_ray(lmgr, p2, p_end, tweight, weight0, color, line_length, cindex, ltype, depth+1, rgen, accum_map, bcube);} // transmitted
			}
			weight *= rweight; // reflected weight
		}
		weight *= (DIFFUSE_REFL*(1.0f - specular) + SPEC_REFL*specular);
	}
	//if (rgen.rand_float() < fabs(weight)/last_weight) return; weight = last_weight; // end the ray (Russian roulette)
	if (fabs(weight) < WEIGHT_THRESH*weight0) return;

	// create reflected ray and make recursive call(s)
	unsigned const num_splits(((depth == 0) ? INIT_RAY_SPLITS : NUM_RAY_SPLITS)[clamp_ltype_range(ltype)]);
	vector3d v_new, v_ref(zero_vector);

	for (unsigned n = 0; n < num_splits; ++n) {
		vector3d const rand_dir(rgen.signed_rand_vector().get_norm());

		if (specular > 0.0 && shine > 1.0 && specular >= rgen.rand_float()) { // specular reflection
			if (v_ref == zero_vector) {
				calc_reflection_angle(dir, v_ref, cnorm);
				v_ref.normalize();
			}
			v_new = (v_ref + rand_dir/sqrt(shine)).get_norm(); // Note: not physically correct
			if (dot_product(v_new, cnorm) < 0.0) continue; // rarely happens?
		}
		else { // random diffuse/Lambertian scatter (cosine distribution using normal-offset sphere)
			v_new = (cnorm + rand_dir).get_norm();
			//assert(dot_product(v_new, cnorm) >= 0.0); // too strong - may fail due to FP rounding
		}
		p2 = p1 + v_new*line_length; // ending point: effectively at infinity
		cast_light_ray(lmgr, cpos, p2, weight/num_splits, weight0, color, line_length, cindex, ltype, depth+1, rgen, accum_map, bcube);
	}
}


struct rt_data {
	unsigned ix, num, job_id, checksum;
	int rseed, ltype;
	bool is_thread, verbose, randomized, is_running;
	cube_t update_bcube;
	lmap_manager_t *lmgr;
	cobj_ray_accum_map_t accum_map;

	rt_data(unsigned i=0, unsigned n=0, int s=1, bool t=0, bool v=0, bool r=0, int lt=0, unsigned jid=0)
		: ix(i), num(n), job_id(jid), checksum(0), rseed(s), ltype(lt), is_thread(t), verbose(v), randomized(r), is_running(0), lmgr(nullptr) {update_bcube.set_to_zeros();}

	void pre_run(rand_gen_t &rgen) {
		assert(lmgr);
		assert(num > 0);
		assert(!is_running);
		is_running = 1;
		rgen.set_state(rseed, 1);
	}
	void post_run() {
		assert(is_running); // can this fail due to race conditions? too strong? remove?
		is_running = 0;
	}
};


template<typename T> class thread_manager_t {

	vector<std::thread> threads;
public:
	vector<T> data; // to be filled in by the caller

	bool is_active() const {return (!threads.empty());}

	bool any_threads_running() const {
		for (auto i = data.begin(); i != data.end(); ++i) {if (i->is_running) return 1;}
		return 0;
	}
	void clear() {
		data.clear();
		threads.clear();
	}
	void create(unsigned num_threads) {
		assert(!is_active());
		data.resize(num_threads);
		threads.resize(num_threads);
	}
	void run(void (*func)(rt_data *)) {
		assert(threads.size() == data.size());
		for (unsigned t = 0; t < threads.size(); ++t) {threads[t] = std::thread(func, (rt_data *)(&data[t]));}
	}
	void join() {
		for (unsigned t = 0; t < threads.size(); ++t) {threads[t].join();}
	}
	void join_and_clear() {join(); clear();}
};

thread_manager_t<rt_data> thread_manager;
lmap_manager_t thread_temp_lmap;

bool indir_lighting_updated() {return (global_lighting_update && (lmap_manager.was_updated || thread_temp_lmap.was_updated));} // only for global updates


void kill_current_raytrace_threads() {

	if (thread_manager.is_active()) { // can't have two running at once, so kill the existing one
		// cancel thread?
		kill_raytrace = 1;
		thread_manager.join_and_clear();
		assert(!thread_manager.is_active());
		kill_raytrace = 0;
	}
}


void update_lmap_from_temp_copy() {

	if (!thread_temp_lmap.was_updated) return; // no updates
	float const blend_weight = 1.0; // TODO: slow blend over time to reduce popping
	lmap_manager.copy_data(thread_temp_lmap, blend_weight);
	thread_temp_lmap.was_updated = 0;
	lmap_manager.was_updated     = 1;
}


void check_for_lighting_finished() { // to be called about once per frame

	if (!thread_manager.is_active()) return; // inactive
	if (thread_manager.any_threads_running()) return; // still running
	thread_manager.join_and_clear(); // clear() or join_and_clear()?
	update_lmap_from_temp_copy();
}


// see https://computing.llnl.gov/tutorials/pthreads/ (for old pthread implementation - now using std::thread)
void launch_threaded_job(unsigned num_threads, void (*start_func)(rt_data *), bool verbose, bool blocking, bool use_temp_lmap, bool randomized, int ltype, unsigned job_id=0) {

	kill_current_raytrace_threads();
	assert(num_threads > 0 && num_threads < 100);
	assert(!keep_beams || num_threads == 1); // could use a mutex instead to make this legal
	bool const single_thread(num_threads == 1);
	if (verbose) {cout << "Computing lighting on " << num_threads << " threads." << endl;}
	thread_manager.create(num_threads);
	vector<rt_data> &data(thread_manager.data);
	if (use_temp_lmap) {thread_temp_lmap.init_from(lmap_manager);}

	for (unsigned t = 0; t < data.size(); ++t) {
		// create a custom lmap_manager_t for each thread then merge them together?
		data[t] = rt_data(t, num_threads, 234323*(t+1), !single_thread, (verbose && t == 0), randomized, ltype, job_id);
		data[t].lmgr = (use_temp_lmap ? &thread_temp_lmap : &lmap_manager);
	}
	if (single_thread && blocking) { // threads disabled
		start_func((rt_data *)(&data[0]));
	}
	else {
		thread_manager.run(start_func);
		if (blocking) {thread_manager.join();}
	}
	if (blocking) {
		if (enable_platform_lights(ltype)) {
			merged_accum_map.clear();
			for (auto i = data.begin(); i != data.end(); ++i) {merged_accum_map.merge(i->accum_map);}
			if (!merged_accum_map.empty()) {merged_accum_map.stats();}
		}
		if (ltype == LIGHTING_COBJ_ACCUM) {
			for (auto i = data.begin(); i != data.end(); ++i) {
				lmap_manager.update_bcube.assign_or_union_with_cube(i->update_bcube); // merge update bounding cubes
			}
		}
		thread_manager.clear();
	}
	//cout << "total rays: " << tot_rays << ", hits: " << num_hits << ", cells touched: " << cells_touched << endl;
	//tot_rays = num_hits = cells_touched = 0;
	// if non-blocking, threads are finished when lmap_manager.was_updated is not set to true during the current frame
}


bool cube_light_src_vect::ray_intersects_any(point const &start_pt, point const &end_pt) const {

	// skip any lines that cross a cube light because otherwise we would doubly add the contribution from that pos/dir
	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->intensity != 0.0 && check_line_clip(start_pt, end_pt, i->bounds.d)) return 1;
	}
	return 0;
}


void trace_one_global_ray(lmap_manager_t *lmgr, point const &pos, point const &pt, colorRGBA const &color, float ray_wt,
	int ltype, bool is_scene_cube, rand_gen_t &rgen, cobj_ray_accum_map_t *accum_map, float line_length)
{
	point const end_pt(pt + (pt - pos).get_norm()*line_length);
	if (is_scene_cube && global_cube_lights.ray_intersects_any(pt, end_pt)) return; // don't double count
	cast_light_ray(lmgr, pos, end_pt, ray_wt, ray_wt, color, line_length, -1, ltype, 0, rgen, accum_map);
}


void trace_ray_block_global_cube(lmap_manager_t *lmgr, cube_t const &bnds, point const &pos, colorRGBA const &color, float ray_wt,
	unsigned nrays, int ltype, unsigned disabled_edges, bool is_scene_cube, bool verbose, bool randomized, rand_gen_t &rgen, cobj_ray_accum_map_t *accum_map)
{
	float const line_length(2.0*get_scene_radius());
	vector3d const ldir((bnds.get_cube_center() - pos).get_norm());
	float proj_area[3] = {0}, tot_area(0.0);

	for (unsigned i = 0; i < 3; ++i) { // adjust the number or weight of rays based on sun/moon position, or simply modify color scale?
		if (disabled_edges & EFLAGS[i][(ldir[i] < 0.0)]) continue; // should this be here, or should we just skip them later?
		unsigned const d0((i+1)%3), d1((i+2)%3);
		vector3d norm(0.0, 0.0, 0.0);
		norm[i]      = 1.0;
		proj_area[i] = fabs((bnds.d[d0][1] - bnds.d[d0][0])*(bnds.d[d1][1] - bnds.d[d1][0])*dot_product(ldir, norm));
		tot_area    += proj_area[i];
	}
	for (unsigned i = 0; i < 3; ++i) {
		if (proj_area[i] == 0.0) continue;
		assert(tot_area > 0.0);
		bool const dir(ldir[i] < 0.0);
		unsigned const d0((i+1)%3), d1((i+2)%3);
		unsigned const num_rays(unsigned(nrays*proj_area[i]/tot_area + 0.5));
		point pt;
		pt[i] = bnds.d[i][dir];
		if (verbose) {cout << "Dim " << i+1 << " of 3, num (this thread): " << num_rays << ", progress (of " << 1+num_rays/1000 << "): 0";}

		if (randomized) {
			for (unsigned s = 0; s < num_rays; ++s) {
				if (kill_raytrace) break;
				if (verbose && ((s%1000) == 0)) {increment_printed_number(s/1000);}
				pt[d0] = rgen.rand_uniform(bnds.d[d0][0], bnds.d[d0][1]);
				pt[d1] = rgen.rand_uniform(bnds.d[d1][0], bnds.d[d1][1]);
				trace_one_global_ray(lmgr, pos, pt, color, ray_wt, ltype, is_scene_cube, rgen, accum_map, line_length);
			}
		}
		else {
			float const len0(bnds.d[d0][1] - bnds.d[d0][0]), len1(bnds.d[d1][1] - bnds.d[d1][0]);
			unsigned const n0(max(1U, unsigned(sqrt((float)num_rays)*len0/len1)));
			unsigned const n1(max(1U, unsigned(sqrt((float)num_rays)*len1/len0)));
			unsigned num(0);

			for (unsigned s0 = 0; s0 < n0; ++s0) {
				if (kill_raytrace) break;
				pt[d0] = bnds.d[d0][0] + (s0 + rgen.rand_uniform(0.0, 1.0))*len0/n0;

				for (unsigned s1 = 0; s1 < n1; ++s1, ++num) {
					if (kill_raytrace) break;
					if (verbose && ((num%1000) == 0)) increment_printed_number(num/1000);
					pt[d1] = bnds.d[d1][0] + (s1 + rgen.rand_uniform(0.0, 1.0))*len1/n1;
					trace_one_global_ray(lmgr, pos, pt, color, ray_wt, ltype, is_scene_cube, rgen, accum_map, line_length);
				}
			}
		}
		if (verbose) {cout << endl;}
	} // for i
}


void trace_ray_block_global_light(rt_data *data, point const &pos, colorRGBA const &color, float weight) {

	if (pos.z < 0.0 || weight == 0.0 || color.alpha == 0.0) return; // below the horizon or zero weight, skip it
	assert(data);
	rand_gen_t rgen;
	data->pre_run(rgen);
	unsigned long long cube_start_rays(0);

	if (GLOBAL_RAYS > 0) {
		float const ray_wt(RAY_WEIGHT*weight*color.alpha/GLOBAL_RAYS);
		assert(ray_wt > 0.0);
		cube_t const bnds(get_scene_bounds());
		trace_ray_block_global_cube(data->lmgr, bnds, pos, color, ray_wt, max(1U, GLOBAL_RAYS/data->num), LIGHTING_GLOBAL, 0, 1, data->verbose, data->randomized, rgen, &data->accum_map);
	}
	for (cube_light_src_vect::const_iterator i = global_cube_lights.begin(); i != global_cube_lights.end(); ++i) {
		if (data->num == 0 || i->num_rays == 0) continue; // disabled
		if (data->verbose) {cout << "Cube volume light source " << (i - global_cube_lights.begin()) << " of " << global_cube_lights.size() << endl;}
		unsigned const num_rays(i->num_rays/data->num);
		float const cube_weight(RAY_WEIGHT*weight*i->intensity/i->num_rays);
		trace_ray_block_global_cube(data->lmgr, i->bounds, pos, color, cube_weight, num_rays, LIGHTING_GLOBAL, i->disabled_edges, 0, data->verbose, data->randomized, rgen, &data->accum_map);
		cube_start_rays += num_rays;
	}
	if (data->verbose) {
		cout << "start rays: " << GLOBAL_RAYS << ", cube_start_rays: " << cube_start_rays << ", total rays: "
			 << tot_rays << ", hits: " << num_hits << ", cells touched: " << cells_touched << endl;
	}
	data->post_run();
}


void trace_ray_block_global(rt_data *data) {

	if (GLOBAL_RAYS == 0 && global_cube_lights.empty()) return; // nothing to do
	// Note: The light color here is white because it will be multiplied by the ambient color later,
	//       and the moon color is generally similar to the sun color so they can be approximated as equal
	float const lfn(CLIP_TO_01(1.0f - 5.0f*(light_factor - 0.4f)));
	if (light_factor >= 0.4 && !combined_gu) trace_ray_block_global_light(data, sun_pos,  WHITE, 1.0-lfn);
	if (light_factor <= 0.6) trace_ray_block_global_light(data, moon_pos, WHITE, lfn);
}


float get_sky_light_ray_weight() {return RAY_WEIGHT/(((float)NPTS)*NRAYS);}

void trace_ray_block_sky(rt_data *data) {

	assert(data);
	rand_gen_t rgen;
	data->pre_run(rgen);
	float const scene_radius(get_scene_radius()), line_length(2.0*scene_radius);
	unsigned long long start_rays(0), cube_start_rays(0);

	if (NPTS > 0 && NRAYS > 0) {
		float const ray_wt(get_sky_light_ray_weight());
		unsigned const block_npts(max(1U, NPTS/data->num));
		vector<point> pts(block_npts);
		vector<vector3d> dirs(NRAYS);

		for (unsigned p = 0; p < block_npts; ++p) {
			do {
				pts[p] = rgen.signed_rand_vector_spherical(1.0).get_norm()*scene_radius; // start the ray here
			} while (pts[p].z < zbottom); // force above zbottom
		}
		sort(pts.begin(), pts.end());
		if (data->verbose) {cout << "Sky light source progress (of " << block_npts << "): 0";}

		for (unsigned p = 0; p < block_npts; ++p) {
			if (kill_raytrace) break;
			if (data->verbose) {increment_printed_number(p);}
			point const &pt(pts[p]);

			for (unsigned r = 0; r < NRAYS; ++r) {
				point const target_pt(X_SCENE_SIZE*rgen.signed_rand_float(), Y_SCENE_SIZE*rgen.signed_rand_float(), rgen.rand_uniform(czmin, czmax));
				dirs[r] = (target_pt - pt).get_norm();
				//dirs[r].z = -fabs(dirs[r].z); // pointing down
			}
			sort(dirs.begin(), dirs.end());

			for (unsigned r = 0; r < NRAYS; ++r) {
				if (kill_raytrace) break;
				if (dot_product(dirs[r], pt) >= 0.0) continue; // can get here when (-Z_SCENE_SIZE, Z_SCENE_SIZE) does not contain (czmin, czmax)
				point const end_pt(pt + dirs[r]*line_length);
				if (sky_cube_lights.ray_intersects_any(pt, end_pt)) continue; // don't double count
				cast_light_ray(data->lmgr, pt, end_pt, ray_wt, ray_wt, WHITE, line_length, -1, LIGHTING_SKY, 0, rgen, &data->accum_map);
				++start_rays;
			}
		}
		if (data->verbose) {cout << endl;}
	}
	for (cube_light_src_vect::const_iterator i = sky_cube_lights.begin(); i != sky_cube_lights.end(); ++i) {
		if (kill_raytrace) break;
		if (data->num == 0 || i->num_rays == 0) continue; // disabled
		unsigned const num_rays(i->num_rays/data->num);
		float const cube_weight(RAY_WEIGHT*i->intensity/i->num_rays);
		if (data->verbose) {cout << "Cube volume light source " << (i - sky_cube_lights.begin()) << " of " << sky_cube_lights.size() << ", progress (of " << 1+num_rays/1000 << "): 0";}
		cube_start_rays += num_rays;

		for (unsigned p = 0; p < num_rays; ++p) {
			if (kill_raytrace) break;
			if (data->verbose && ((p%1000) == 0)) {increment_printed_number(p/1000);}
			point const pt(rgen.gen_rand_cube_point(i->bounds));
			vector3d dir(rgen.signed_rand_vector_spherical().get_norm()); // need high quality distribution
			dir.z = -fabs(dir.z); // make sure z is negative since this is supposed to be light from the sky
			point const end_pt(pt + dir*line_length);
			cast_light_ray(data->lmgr, pt, end_pt, cube_weight, cube_weight, i->color, line_length, -1, LIGHTING_SKY, 0, rgen, &data->accum_map);
		}
		if (data->verbose) {cout << endl;}
	}
	if (data->verbose) {
		cout << "start rays: " << start_rays << ", cube start rays: " << cube_start_rays << ", total rays: " << tot_rays
			 << ", hits: " << num_hits << ", cells touched: " << cells_touched << endl;
	}
	data->checksum = rgen.rand();
	data->post_run();
}


bool is_correct_accum_cobj(coll_obj const &cobj, cube_t const bcube) {
	return (cobj.is_update_light_platform() && bcube.intersects(cobj.get_platform_max_bcube()));
}
coll_obj &find_accum_cobj(unsigned id, cobj_ray_accum_t const &arc) {

	cube_t bcube(arc.get_bcube(1));
	coll_obj &cobj(coll_objects.get_cobj(id));
	if (is_correct_accum_cobj(cobj, bcube)) return cobj;

	for (cobj_id_set_t::const_iterator i = coll_objects.platform_ids.begin(); i != coll_objects.platform_ids.end(); ++i) {
		coll_obj &cobj(coll_objects.get_cobj(*i));
		if (is_correct_accum_cobj(cobj, bcube)) return cobj;
	}
	cerr << "Error: Failed to find platform cobj with ID " << id << endl;
	assert(0);
	return cobj; // never gets here
}

void trace_ray_block_cobj_accum(rt_data *data) {

	assert(data);
	rand_gen_t rgen;
	data->pre_run(rgen);
	float const line_length(2.0*get_scene_radius()), ray_wt(get_sky_light_ray_weight()); // Note: weight assumes not using cube sky lights

	for (auto i = merged_accum_map.begin(); i != merged_accum_map.end(); ++i) {
		coll_obj &cobj(find_accum_cobj(i->first, i->second));
		cobj.unexpand_from_platform_max_bounds(); // unexpand if it was expanded

		// round robin distribute rays across threads
		for (auto r = (i->second.rays.begin() + data->ix); r < i->second.rays.end(); r += data->num) {
			if (kill_raytrace) break; // not needed?
			assert(r->weight > 0.0);
			float const weight0(ray_wt ? ray_wt : r->weight);
			cast_light_ray(data->lmgr, r->pos, r->get_p2(line_length), r->weight, weight0, r->get_color(), line_length, -1, LIGHTING_COBJ_ACCUM, 0, rgen, nullptr, nullptr);
		}
	}
	data->post_run();
}

void trace_ray_block_cobj_accum_single_update(rt_data *data) {

	assert(data);
	rand_gen_t rgen;
	data->pre_run(rgen);
	unsigned const cid(data->job_id);
	coll_obj &cobj(coll_objects.get_cobj(cid));
	assert(cobj.is_update_light_platform());
	vector3d const platform_delta(platforms.get_cobj_platform(cobj).get_last_delta());
	if (platform_delta == zero_vector) return; // platform hasn't moved - why are we here?
	data->update_bcube.set_to_zeros();
	cube_t const prev_bcube(cobj - platform_delta); // previous frame's position of cobj
	float const line_length(2.0*get_scene_radius()), ray_wt(get_sky_light_ray_weight()); // Note: weight assumes not using cube sky lights
	auto it(merged_accum_map.find(cid));
	
	if (it == merged_accum_map.end()) { // not the correct cobj, search for correct accum map entry
		for (auto i = merged_accum_map.begin(); i != merged_accum_map.end(); ++i) {
			if (is_correct_accum_cobj(cobj, i->second.get_bcube(1))) {it = i; break;}
		}
		assert(it != merged_accum_map.end());
	}
	// round robin distribute rays across threads
	for (auto r = (it->second.rays.begin() + data->ix); r < it->second.rays.end(); r += data->num) {
		assert(r->weight > 0.0);
		point const end_pt(r->get_p2(line_length));
		bool const cur_hit(check_line_clip(r->pos, end_pt, cobj.d)), prev_hit(check_line_clip(r->pos, end_pt, prev_bcube.d));
		if (cur_hit == prev_hit) continue; // no change in hit status
		float const weight(r->weight*(cur_hit ? -1.0 : 1.0)); // if ray is newly blocked, subtract its contribution by negating its weight
		// Note: cobj is ignored here because it can't be in both the prev and cur position at the same time, and temporarily moving it isn't thread safe
		cast_light_ray(data->lmgr, r->pos, end_pt, weight, (ray_wt ? ray_wt : r->weight), r->get_color(), line_length, cid, LIGHTING_COBJ_ACCUM, 0, rgen, nullptr, &data->update_bcube);
	}
	data->post_run();
}


void ray_trace_local_light_source(lmap_manager_t *lmgr, light_source const &ls, float line_length, unsigned num_rays, rand_gen_t &rgen, int ltype, unsigned N_RAYS) {

	colorRGBA lcolor(ls.get_color());
	if (N_RAYS == 0 || lcolor.alpha == 0.0) return; // nothing to do
	bool const line_light(ls.is_line_light());
	point const &lpos(ls.get_pos()), &lpos2(ls.get_pos2());
	vector3d delta; // only nonzero for line lights
	if (line_light && num_rays > 1) {delta = (lpos2 - lpos)/float(num_rays-1);}
	float const radius(ls.get_radius()), r_inner(ls.get_r_inner()), ray_wt(1000.0*lcolor.alpha*radius/N_RAYS);
	assert(ray_wt > 0.0);

	if (ls.get_is_cube_light()) {
		assert(!line_light);
		cube_t const cube(lpos, lpos2);
		float tot_area(0.0), side_area[3];
		
		for (unsigned dim = 0; dim < 3; ++dim) {
			unsigned const d1((dim+1)%3), d2((dim+2)%3);
			side_area[dim] = (fabs(cube.d[d1][1] - cube.d[d1][0])*fabs(cube.d[d2][1] - cube.d[d2][0]));

			for (unsigned dir = 0; dir < 2; ++dir) {
				if (!(ls.get_cube_eflags() & EFLAGS[dim][dir])) {tot_area += side_area[dim];}
			}
		}
		assert(tot_area > 0.0);
		//cout << TXT(tot_area) << TXT(radius) << TXT(ray_wt) << TXT(num_rays) << TXT(N_RAYS) << endl;

		for (unsigned dim = 0; dim < 3; ++dim) {
			unsigned const d1((dim+1)%3), d2((dim+2)%3);
			unsigned const side_rays(unsigned(num_rays*side_area[dim]/tot_area)); // number of rays for each of the two sides
			if (side_rays == 0) continue;

			for (unsigned dir = 0; dir < 2; ++dir) {
				if (ls.get_cube_eflags() & EFLAGS[dim][dir]) continue; // this surface is disabled
				vector3d normal(zero_vector);
				normal[dim] = (dir ? 1.0 : -1.0);
				point start_pt;
				start_pt[dim] = cube.d[dim][dir] + 1.0E-5*radius*normal[dim]; // move slightly away from cube edge
				//cout << TXT(dim) << TXT(dir) << TXT(d1) << TXT(d2) << TXT(side_area[dim]) << TXT(side_rays) << endl;

				for (unsigned n = 0; n < side_rays; ++n) {
					if (kill_raytrace) break;
					vector3d dir(rgen.signed_rand_vector_spherical(1.0).get_norm());
					if (dot_product(dir, normal) < 0.0) {dir.negate();}
					start_pt[d1] = rgen.rand_uniform(cube.d[d1][0], cube.d[d1][1]);
					start_pt[d2] = rgen.rand_uniform(cube.d[d2][0], cube.d[d2][1]);
					point const end_pt(start_pt + dir*line_length);
					cast_light_ray(lmgr, start_pt, end_pt, ray_wt, ray_wt, lcolor, line_length, -1, ltype, 0, rgen, nullptr); // init_cobj not used here
				} // for n
			} // for dir
		} // for dim
		return;
	} // end cube light case

	int init_cobj(-1);
	check_coll_line(lpos, lpos2, init_cobj, -1, 1, 2); // find most opaque (max alpha) containing object
	assert(init_cobj < (int)coll_objects.size());

	for (unsigned n = 0; n < num_rays; ++n) {
		if (kill_raytrace) break;
		vector3d dir;
		float weight(0.0);

		for (unsigned tries = 0; tries < 10; ++tries) { // make up to 10 attempts
			dir    = rgen.signed_rand_vector_spherical(1.0).get_norm();
			weight = ray_wt*ls.get_dir_intensity(-dir);
			if (weight > 0.0) break; // success
		}
		if (weight == 0.0) continue; // failed to choose a dir
		point start_pt;

		if (init_cobj >= 0 && r_inner == 0.0) { // use a volume light source
			coll_obj cobj(coll_objects[init_cobj]);
			
			while (1) { // generate a radom point within cobj's bounding cube
				start_pt = rgen.gen_rand_cube_point(cobj);
				if (cobj.line_intersect(start_pt, start_pt)) break; // intersection (point contained)
			}
		}
		else { // use a point light source
			start_pt = lpos;

			if (r_inner > 0.0) {
				// move r_inner away from the light source
				// necessary for sources contained in more than one cobj (like the lamps in mapx)
				// move in a different dir to simulate emission along the surface and avoid single point radial grid artifacts
				vector3d const move_dir(rgen.signed_rand_vector_spherical(1.0).get_norm());
				bool const invert(dot_product(dir, move_dir) < 0);
				start_pt += move_dir*(invert ? -r_inner : r_inner);
			}
			if (line_light) {start_pt += n*delta;} // fixed spacing along the length of the line
		}
		point const end_pt(start_pt + dir*line_length);
		cast_light_ray(lmgr, start_pt, end_pt, weight, weight, lcolor, line_length, init_cobj, ltype, 0, rgen, nullptr);
	} // for n
}


void trace_ray_block_local(rt_data *data) {

	assert(data);
	if (LOCAL_RAYS == 0) return; // nothing to do
	rand_gen_t rgen;
	data->pre_run(rgen);
	float const line_length(2.0*get_scene_radius());
	
	if (data->verbose) {
		cout << "Local light sources progress (of " << light_sources_a.size() << "): 0";
		cout.flush();
	}
	for (unsigned i = 0; i < light_sources_a.size(); ++i) {
		if (data->verbose) {increment_printed_number(i);}
		unsigned const light_nrays(light_sources_a[i].get_num_rays()), NRAYS(light_nrays ? light_nrays : LOCAL_RAYS), num_rays(max(1U, NRAYS/data->num));
		ray_trace_local_light_source(data->lmgr, light_sources_a[i], line_length, num_rays, rgen, data->ltype, NRAYS);
	}
	if (data->verbose) {cout << endl;}
	data->post_run();
}


void trace_ray_block_dynamic(rt_data *data) {

	assert(data);
	if (DYNAMIC_RAYS == 0) return; // nothing to do
	light_volume_local const &lvol(get_local_light_volume(data->ltype));
	vector<unsigned> const &dlight_ixs(indir_dlight_group_manager.get_dlight_ixs_for_tag_ix(lvol.get_tag_ix()));
	assert(!dlight_ixs.empty());
	rand_gen_t rgen;
	data->pre_run(rgen);
	float const max_line_length(2.0*get_scene_radius());
	
	for (auto i = dlight_ixs.begin(); i != dlight_ixs.end(); ++i) {
		assert(*i < light_sources_d.size());
		light_source_trig const &ls(light_sources_d[*i]);
		//if (!ls.is_enabled()) continue; // error?
		float const line_length(min(4.0f*ls.get_radius(), max_line_length)); // limit ray length to improve perf
		unsigned const light_nrays(ls.get_num_rays()), NRAYS(light_nrays ? light_nrays : DYNAMIC_RAYS), num_rays(max(1U, NRAYS/data->num));
		ray_trace_local_light_source(nullptr, ls, line_length, num_rays, rgen, data->ltype, NRAYS); // lmgr is unused, so leave it as null
	}
	data->post_run();
}


typedef void (*ray_trace_func)(rt_data *);
ray_trace_func const rt_funcs[NUM_LIGHTING_TYPES] = {trace_ray_block_sky, trace_ray_block_global, trace_ray_block_local, trace_ray_block_cobj_accum, trace_ray_block_dynamic};


void compute_ray_trace_lighting(unsigned ltype, bool verbose) {

	bool const dynamic(is_ltype_dynamic(ltype));
	unsigned const c_ltype(clamp_ltype_range(ltype));
	assert(c_ltype < NUM_LIGHTING_TYPES);
	const char *fn(lighting_file[c_ltype]);

	if (!dynamic && read_light_files[c_ltype]) {
		if (c_ltype == LIGHTING_COBJ_ACCUM) {
			merged_accum_map.open_and_read(fn, 0);

			if (store_cobj_accum_lighting_as_blocked) {
				timer_t t("Cobj Accum Lighting");
				launch_threaded_job(NUM_THREADS, rt_funcs[c_ltype], verbose, 1, 0, 0, ltype); // update fully blocked lighting with currently blocked portion
			}
		}
		else {lmap_manager.read_data_from_file(fn, c_ltype);}
	}
	else {
		if (c_ltype != LIGHTING_LOCAL && !dynamic) {cout << X_SCENE_SIZE << " " << Y_SCENE_SIZE << " " << Z_SCENE_SIZE << " " << czmin << " " << czmax << endl;}
		all_models.build_cobj_trees(1);
		if (enable_platform_lights(ltype)) {pre_rt_bvh_build_hook();}
		launch_threaded_job(NUM_THREADS, rt_funcs[c_ltype], verbose, 1, 0, 0, ltype);
		if (enable_platform_lights(ltype)) {post_rt_bvh_build_hook();}
	}
	if (!dynamic && write_light_files[c_ltype]) {
		if (c_ltype == LIGHTING_COBJ_ACCUM) {
			merged_accum_map.open_and_write(fn, 0);
			// if writing both the cobj accum file and the sky lighting file, and not storing sky lighting as blocked,
			// we need to rewrite sky lighting to include unblocked cobj accum light
			if (!store_cobj_accum_lighting_as_blocked && write_light_files[LIGHTING_SKY]) {
				lmap_manager.write_data_to_file(lighting_file[LIGHTING_SKY], LIGHTING_SKY);
			}
		}
		else {lmap_manager.write_data_to_file(fn, c_ltype);}
	}
}


bool pre_lighting_update() {
	if (!lmap_manager.is_allocated()) return 0; // too early
	tot_rays = num_hits = cells_touched = 0;
	all_models.build_cobj_trees(0);
	return 1;
}

void check_update_global_lighting(unsigned lights) {

	if (!global_lighting_update || !(read_light_files[LIGHTING_GLOBAL] || write_light_files[LIGHTING_GLOBAL])) return;
	if (!(lights & (SUN_SHADOW | MOON_SHADOW))) return;
	if (GLOBAL_RAYS == 0 && global_cube_lights.empty()) return; // nothing to do
	if (!pre_lighting_update()) return; // lmap is not yet allocated
	// Note: we could check if the sun/moon is visible, but it might have been visible previously and now is not, and in that case we still need to update lighting
	no_stat_moving = 1; // disable static moving cobjs for async updates, which aren't thread safe because the BVH is rebuilt every frame; no need to set back after first frame
	lmap_manager.clear_lighting_values(LIGHTING_GLOBAL);
	bool const reserve_thread((int)NUM_THREADS >= omp_get_max_threads()); // reserve a thread for rendering if needed
	launch_threaded_job(max(1U, NUM_THREADS-reserve_thread), rt_funcs[LIGHTING_GLOBAL], 0, 0, lighting_update_offline, 0, LIGHTING_GLOBAL);
}

void check_all_platform_cobj_lighting_update() {

	if (merged_accum_map.empty()) return; // updates not enabled
	if (!pre_lighting_update())   return; // lmap is not yet allocated
	bool const prev_was_updated(lmap_manager.was_updated);
	lmap_manager.was_updated = 0; // clear and check if it gets set again
	cube_t &lm_bc(lmap_manager.update_bcube);

	for (cobj_id_set_t::const_iterator i = coll_objects.platform_ids.begin(); i != coll_objects.platform_ids.end(); ++i) {
		coll_obj &cobj(coll_objects.get_cobj(*i));
		if (platforms.get_cobj_platform(cobj).get_last_delta() == zero_vector) continue; // not moving
		if (!cobj.is_update_light_platform()) continue; // no updates
		launch_threaded_job(NUM_THREADS, trace_ray_block_cobj_accum_single_update, 0, 1, 0, 0, LIGHTING_COBJ_ACCUM, *i); // blocking, on all threads, using cobj_id as job_id
	}
	if (lmap_manager.was_updated && !lm_bc.is_zero_area()) {
		lmap_manager.was_updated = 0; // unset to enable multi-threaded updates (though it doesn't seem to matter much)
		int const x1(max(get_xpos_round_down(lm_bc.d[0][0]), 0)), x2(min(get_ypos_round_down(lm_bc.d[0][1])+1, MESH_X_SIZE));
		int const y1(max(get_xpos_round_down(lm_bc.d[1][0]), 0)), y2(min(get_ypos_round_down(lm_bc.d[1][1])+1, MESH_Y_SIZE));
		int const z1(max(get_zpos(lm_bc.d[2][0]), 0)), z2(min(get_zpos(lm_bc.d[2][1])+1, MESH_SIZE[2]));
		if (x1 < x2 && y1 < y2 && z1 < z2) {update_smoke_indir_tex_range(x1, x2, y1, y2, z1, z2);}
		lm_bc.set_to_zeros(); // clear
	}
	lmap_manager.was_updated = prev_was_updated; // restore previous value
}


// lmap_manager_t


bool lmap_manager_t::read_data_from_file(char const *const fn, int ltype) {

	assert(fn != nullptr);
	binary_file_reader reader;
	if (!reader.open(fn)) return 0;
	cout << "Reading lighting file from " << fn << endl;
	unsigned data_size(0);
	if (!reader.read(&data_size, sizeof(unsigned), 1)) return 0;

	if (data_size != vldata_alloc.size()) {
		cerr << "Error: Lighting file " << fn << " data size of " << data_size
			 << " does not equal the expected size of " << vldata_alloc.size() << ". Ignoring file." << endl;
		return 0;
	}
	unsigned const sz = lmcell::get_dsz(ltype);
	vector<float> data(data_size*sz);
	unsigned pos(0);

	if (!reader.read(&data.front(), sizeof(float), data.size())) {
		cerr << "Error reading data from ligthing file " << fn << endl;
		return 0;
	}
	for (lmcell &c : vldata_alloc) {
		float *ptr(c.get_offset(ltype));
		for (unsigned n = 0; n < sz; ++n) {ptr[n] = data[pos++];}
	}
	assert(pos == data.size());
	return 1;
}


bool lmap_manager_t::write_data_to_file(char const *const fn, int ltype) const {

	if (fn == nullptr || strcmp(fn, "''") == 0 || strcmp(fn, "\"\"") == 0) return 0; // don't write
	binary_file_writer writer;
	if (!writer.open(fn)) return 0;
	cout << "Writing lighting file to " << fn << endl;
	unsigned const data_size((unsigned)vldata_alloc.size()); // should be size_t?
	if (!writer.write(&data_size, sizeof(unsigned), 1)) return 0;
	unsigned const sz(lmcell::get_dsz(ltype));

	for (vector<lmcell>::const_iterator i = vldata_alloc.begin(); i != vldata_alloc.end(); ++i) { // const_iterator?
		if (!writer.write(i->get_offset(ltype), sizeof(float), sz)) {
			cerr << "Error writing data to ligthing file " << fn << endl;
			return 0;
		}
	}
	return 1;
}


void lmap_manager_t::clear_lighting_values(int ltype) {

	assert(ltype < NUM_LIGHTING_TYPES && !is_ltype_dynamic(ltype));
	unsigned const num(lmcell::get_dsz(ltype));

	for (vector<lmcell>::iterator i = vldata_alloc.begin(); i != vldata_alloc.end(); ++i) {
		float *color(i->get_offset(ltype));
		for (unsigned j = 0; j < num; ++j) {color[j] = 0.0;}
	}
}

