// 3D World - Ray Tracing Code
// by Frank Gennari
// 2/14/10
#include <pthread.h>
#include "3DWorld.h"
#include "lightmap.h"
#include "gameplay.h"


float const RAY_WEIGHT    = 4.0E5;
float const WEIGHT_THRESH = 0.01;
float const DIFFUSE_REFL  = 0.9; // 90%  diffuse  reflectivity
float const SPEC_REFL     = 1.0; // 100% specular reflectivity
float const SNOW_ALBEDO   = 0.9;
float const ICE_ALBEDO    = 0.8;

bool keep_lasers(0); // debugging mode
unsigned NPTS(50000), NRAYS(40000), LOCAL_RAYS(1000000), GLOBAL_RAYS(1000000), NUM_THREADS(1);
unsigned long long tot_rays(0), num_hits(0), cells_touched(0);
unsigned const NUM_RAY_SPLITS [NUM_LIGHTING_TYPES] = {1, 1, 1}; // sky, global, local
unsigned const INIT_RAY_SPLITS[NUM_LIGHTING_TYPES] = {1, 4, 1}; // sky, global, local

extern bool has_snow;
extern int read_light_files[], write_light_files[], display_mode, DISABLE_WATER;
extern float light_int_scale[], ztop, water_plane_z, temperature, snow_depth, indir_light_exp, first_ray_weight;
extern char *lighting_file[];
extern point sun_pos, moon_pos;
extern vector<light_source> light_sources;
extern vector<coll_obj> coll_objects;
extern vector<laser_beam> lasers;
extern lmap_manager_t lmap_manager;
extern vector<cube_light_source> global_cube_lights;


float get_scene_radius() {return sqrt(2.0*(X_SCENE_SIZE*X_SCENE_SIZE + Y_SCENE_SIZE*Y_SCENE_SIZE + Z_SCENE_SIZE*Z_SCENE_SIZE));}
float get_step_size()    {return 0.3*(DX_VAL + DY_VAL + DZ_VAL);}


void increment_printed_number(unsigned num) {

	for (unsigned n = max(num, 1U); n > 0; n /= 10) cout << "\b";
	cout << (num+1);
	cout.flush();
}


void add_path_to_lmcs(point p1, point p2, float weight, colorRGBA const &color, int ltype, bool first_pt) {

	if (first_pt && ltype == LIGHTING_GLOBAL) weight *= first_ray_weight; // lower weight - handled by direct illumination
	colorRGBA const cw(color*weight);
	float const dist(p2p_dist(p1, p2)); // dist can be 0
	unsigned const nsteps(1 + unsigned(dist/get_step_size())); // round up
	vector3d const step((p2 - p1)/nsteps); // at least two points
	if (!first_pt) p1 += step; // move past the first step so we don't double count

	for (unsigned s = 0; s < nsteps+first_pt; ++s) {
		lmcell *lmc(lmap_manager.get_lmcell(p1));
		
		if (lmc != NULL) { // could use a pthread_mutex_t here, but it seems too slow
			float *color(lmc->get_offset(ltype));
			ADD_LIGHT_CONTRIB(cw, color);
			if (ltype != LIGHTING_LOCAL) color[3] += weight;
		}
		p1 += step;
	}
	cells_touched += nsteps;
	++num_hits;
}


void cast_light_ray(point p1, point p2, float weight, float weight0, colorRGBA color,
					float line_length, int ignore_cobj, int ltype, bool first_pt)
{
	assert(!is_nan(p1) && !is_nan(p2));
	++tot_rays;

	// find intersection point with scene cobjs
	point orig_p1(p1);
	if (!do_line_clip_scene(p1, p2, min(zbottom, czmin), max(ztop, czmax)) || is_under_mesh(p1)) return;
	int cindex(-1), xpos, ypos;
	point cpos(p2);
	vector3d cnorm;
	float t, zval;
	bool snow_coll(0), ice_coll(0), water_coll(0), mesh_coll(0);
	vector3d const dir((p2 - p1).get_norm());
	bool coll(check_coll_line_exact(p1, p2, cpos, cnorm, cindex, 0.0, ignore_cobj, 1, 0, 1)); // fast=1
	assert(coll ? (cindex >= 0 && cindex < (int)coll_objects.size()) : (cindex == -1));

	// p1 starts inside a cobj, so find the intersection point with a reverse ray
	if (first_pt && coll && cpos == orig_p1) {
		if (coll_objects[cindex].line_int_exact(p2, p1, t, cnorm, 0.0, 1.0)) {
			cpos = (p2 + (p1 - p2)*t);
			cast_light_ray(cpos, p2, weight, weight0, color, line_length, cindex, ltype, 0);
			return;
		}
	}

	// find intersection point with mesh (approximate)
	if ((display_mode & 0x01) && !coll && p1.z != p2.z && line_intersect_mesh(p1, p2, xpos, ypos, zval, 0, 0)) {
		assert(!point_outside_mesh(xpos, ypos));
		
		if (!is_mesh_disabled(xpos, ypos)) {
			if (p2.z > p1.z) return; // starts under mesh = bad
			cpos  = (p1 + (p2 - p1)*(zval + SMALL_NUMBER - p1.z)/(p2.z - p1.z));
			cnorm = vertex_normals[ypos][xpos];
			coll  = 1;
			mesh_coll = 1;
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
	if ( coll) p2 = cpos;
	if (keep_lasers && p1 != p2) lasers.push_back(laser_beam(!coll, 1, p1, p2, color, 0.1*weight)); // testing
	if (!coll) return; // more efficient to do this up here and let a reverse ray from the sky light this path

	// walk from p1 to p2, adding light to all lightmap cells encountered
	add_path_to_lmcs(p1, p2, weight, color, ltype, first_pt);
	//if (!coll)    return;
	if (p1 == p2) return; // line must have started inside a cobj - this is bad, but what can we do?

	// update ray weight based on object material properties
	float specular(0.0);

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
		float const depth(water_plane_z - cpos.z);
		float const dist(-2*depth*delta.mag()/delta.z); // multiply by 2 to account for both directions
		assert(dist >= 0.0);
		colorRGBA water_color(WATER_C);
		water_color.alpha = 1.0;  // make solid for volume attenuation (not surface)
		water_color      *= 0.95; // refracted component
		atten_by_water_depth(&(water_color.red), 0.8*dist); // assumes dist scaling for outside water
		weight *= water_color.get_luminance();
		color   = color.modulate_with(water_color);
	} // Note: no else
	if (snow_coll) {
		weight  *= SNOW_ALBEDO; // white
		specular = 0.8;
	}
	else if (ice_coll) {
		weight  *= ICE_ALBEDO*ICE_C.get_luminance();
		color    = color.modulate_with(ICE_C);
		specular = 0.9;
	}
	else if (mesh_coll) { // collision with mesh
		colorRGBA const lc(get_landscape_texture_color(xpos, ypos));
		weight *= DIFFUSE_REFL*lc.get_luminance(); // 90% diffuse reflectivity
		color   = color.modulate_with(lc);
	}
	else { // collision with cobj
		assert(cindex >= 0);
		coll_obj const &cobj(coll_objects[cindex]);
		float luminance(1.0), alpha(1.0);
		get_lum_alpha(cobj.cp.color, cobj.cp.tid, luminance, alpha); // use alpha?
		specular = cobj.cp.specular;
		weight  *= luminance;
		color    = color.modulate_with(cobj.cp.color);
		if (cobj.cp.tid >= 0) color = color.modulate_with(texture_color(cobj.cp.tid));

		// calculate refracted ray
		if (alpha < 1.0) { // semi-transparent
			float rweight(alpha);

			if (cobj.cp.refract_ix != 1.0) {
				rweight = get_reflected_weight(get_fresnel_reflection(dir, cnorm*-1, 1.0, cobj.cp.refract_ix), alpha);
			}
			float tweight((1.0 - rweight)*weight); // refracted weight
			
			if (tweight > WEIGHT_THRESH*weight0) {
				bool no_transmit(0);

				if (cobj.cp.refract_ix != 1.0) { // refracted
					vector3d v_refract, v_refract2;
					
					if (calc_refraction_angle(dir, v_refract, cnorm, 1.0, cobj.cp.refract_ix)) {
						point const enter_pt(p2);
						p_end = (p2 + v_refract*line_length);
						vector3d cnorm2;

						// test for collision with reversed ray to get the other intersection point
						if (cobj.line_int_exact(p_end, p2, t, cnorm2, 0.0, 1.0)) { // not sure what to do if fails or tmax >= 1.0
							point const p_int(p_end + (p2 - p_end)*t);
							if (!dist_less_than(p2, p_int, get_step_size())) add_path_to_lmcs(p2, p_int, weight, color, ltype, first_pt);
							
							if (calc_refraction_angle(v_refract, v_refract2, cnorm2*-1, cobj.cp.refract_ix, 1.0)) {
								p2    = p_int;
								p_end = p2 + v_refract2*line_length;
								tweight    *= cobj.get_light_transmit(enter_pt, p_int); // can we use p2p_dist(enter_pt, p_int) directly?
								no_transmit = !(tweight > WEIGHT_THRESH*weight0);
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
				if (!no_transmit) cast_light_ray(p2, p_end, tweight, weight0, color, line_length, cindex, ltype, 0); // transmitted
			}
			weight *= rweight; // reflected weight
		}
		weight *= (DIFFUSE_REFL*(1.0 - specular) + SPEC_REFL*specular);
	}
	if (weight < WEIGHT_THRESH*weight0) return;

	// create reflected ray and make recursive call(s)
	vector3d v_ref;
	calc_reflection_angle(dir, v_ref, cnorm);
	v_ref.normalize();
	unsigned const num_splits(first_pt ? INIT_RAY_SPLITS[ltype] : NUM_RAY_SPLITS[ltype]);
	
	for (unsigned n = 0; n < num_splits; ++n) {
		vector3d new_v(signed_rand_vector().get_norm()); // add random diffuse scatter
		if (dot_product(new_v, cnorm) < 0.0) new_v.negate(); // make in same direction as normal
		vector3d const v_new((v_ref*specular + new_v*(1.0 - specular)).get_norm());
		p2 = p1 + v_new*line_length; // ending point: effectively at infinity
		cast_light_ray(cpos, p2, weight/num_splits, weight0, color, line_length, cindex, ltype, 0);
	}
}


struct rt_data {
	unsigned ix, num;
	int rseed;
	bool is_thread, verbose;

	rt_data(unsigned i=0, unsigned n=0, int s=0, bool t=0, bool v=0) : ix(i), num(n), rseed(s), is_thread(t), verbose(v) {}
};


// see https://computing.llnl.gov/tutorials/pthreads/
void launch_threaded_job(unsigned num_threads, void *(*start_func)(void *), bool verbose) {

	assert(num_threads > 0 && num_threads < 100);
	assert(!keep_lasers || NUM_THREADS == 1); // could use a pthread_mutex_t instead to make this legal
	vector<rt_data> data(num_threads);
	bool const single_thread(num_threads == 1);

	for (unsigned t = 0; t < data.size(); ++t) {
		// create a custom lmap_manager_t for each thread then merge them together?
		data[t] = rt_data(t, num_threads, 234323*t, !single_thread, (verbose && single_thread));
	}
	if (single_thread) { // pthreads disabled
		start_func((void *)(&data[0]));
		return;
	}
	vector<pthread_t> threads(num_threads);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (unsigned t = 0; t < num_threads; ++t) {
		int const rc(pthread_create(&threads[t], &attr, start_func, (void *)(&data[t]))); 
		
		if (rc) {
			cout << "Error: Return code from pthread_create() is " << rc << endl;
			assert(0);
		}
	}
	pthread_attr_destroy(&attr);

	for (unsigned t = 0; t < num_threads; ++t) {
		int const rc(pthread_join(threads[t], NULL));
		
		if (rc) {
			cout << "Error: return code from pthread_join() is " << rc << endl;
			assert(0);
		}
	}
}


void trace_ray_block_global_light(void *ptr, point const &pos, colorRGBA const &color, float weight) {

	if (pos.z < 0.0) return; // below the horizon, skip it
	assert(ptr);
	rt_data *data((rt_data *)ptr);
	cout << "Starting on thread " << data->ix << endl;
	assert(data->num > 0);
	if (data->is_thread) srand(data->rseed);
	float const scene_radius(get_scene_radius()), line_length(2.0*scene_radius);
	float const ray_wt(2.0E5*weight*color.alpha/GLOBAL_RAYS);
	if (ray_wt == 0.0) return; // error?
	unsigned const num_tot_rays(max(1U, GLOBAL_RAYS/data->num));
	cube_t const bnds(-X_SCENE_SIZE, X_SCENE_SIZE, -Y_SCENE_SIZE, Y_SCENE_SIZE, min(zbottom, czmin), max(ztop, czmax));
	vector3d const ldir((bnds.get_cube_center() - pos).get_norm());
	float proj_area[3], tot_area(0.0);

	for (unsigned i = 0; i < 3; ++i) {
		unsigned const d0((i+1)%3), d1((i+2)%3);
		vector3d norm(0.0, 0.0, 0.0);
		norm[i] = 1.0;
		proj_area[i] = fabs((bnds.d[d0][1] - bnds.d[d0][0])*(bnds.d[d1][1] - bnds.d[d1][0])*dot_product(ldir, norm));
		tot_area    += proj_area[i];
	}
	for (unsigned i = 0; i < 3; ++i) {
		unsigned const d0((i+1)%3), d1((i+2)%3);
		unsigned const num_rays(unsigned(num_tot_rays*proj_area[i]/tot_area + 0.5));
		bool const dir(ldir[i] < 0.0);
		float const len0(bnds.d[d0][1] - bnds.d[d0][0]), len1(bnds.d[d1][1] - bnds.d[d1][0]);
		unsigned const n0(max(1U, unsigned(sqrt((float)num_rays)*len0/len1)));
		unsigned const n1(max(1U, unsigned(sqrt((float)num_rays)*len1/len0)));
		point pt;
		pt[i] = bnds.d[i][dir];
		if (data->verbose) cout << "Dim " << i+1 << " of 3" << endl;

		for (unsigned s0 = 0; s0 < n0; ++s0) {
			pt[d0] = bnds.d[d0][0] + s0*len0/n0;

			for (unsigned s1 = 0; s1 < n1; ++s1) {
				pt[d1] = bnds.d[d1][0] + s1*len1/n1;
				point const end_pt(pt + (pt - pos).get_norm()*line_length);
				cast_light_ray(pos, end_pt, ray_wt, ray_wt, color, line_length, -1, LIGHTING_GLOBAL, 1);
			}
		}
	}
	if (data->verbose) {
		cout << "start rays: " << GLOBAL_RAYS << ", total rays: " << tot_rays << ", hits: " << num_hits << ", cells touched: " << cells_touched << endl;
	}
}


void *trace_ray_block_global(void *ptr) {

	if (GLOBAL_RAYS == 0) return 0; // nothing to do
	// FIXME: what about universe mode/combined_gu?
	// Note1: The light color here is white because it will be multiplied by the ambient color later,
	//        and the moon color is generally similar to the sun color so they can be approximated as equal
	float const lfn(CLIP_TO_01(1.0f - 5.0f*(light_factor - 0.4f)));
	if (light_factor >= 0.4) trace_ray_block_global_light(ptr, sun_pos,  WHITE, 1.0-lfn);
	if (light_factor <= 0.6) trace_ray_block_global_light(ptr, moon_pos, WHITE, lfn);
	return 0;
}


void *trace_ray_block_sky(void *ptr) {

	assert(ptr);
	rt_data *data((rt_data *)ptr);
	cout << "Starting on thread " << data->ix << endl;
	assert(data->num > 0);
	if (data->is_thread) srand(data->rseed);
	float const scene_radius(get_scene_radius()), line_length(2.0*scene_radius);
	unsigned long long start_rays(0), cube_start_rays(0);

	if (NPTS > 0 && NRAYS > 0) {
		float const ray_wt(RAY_WEIGHT/(((float)NPTS)*NRAYS));
		unsigned const block_npts(max(1U, NPTS/data->num));
		vector<point> pts(block_npts);
		vector<vector3d> dirs(NRAYS);

		for (unsigned p = 0; p < block_npts; ++p) {
			pts[p] = signed_rand_vector_spherical(1.0, 0).get_norm()*scene_radius; // start the ray here
		}
		sort(pts.begin(), pts.end());
		if (data->verbose) cout << "Sky light source points: " << block_npts << ", progress: 0";

		for (unsigned p = 0; p < block_npts; ++p) {
			if (data->verbose) increment_printed_number(p);
			point const &pt(pts[p]);

			for (unsigned r = 0; r < NRAYS; ++r) {
				point const target_pt(X_SCENE_SIZE*signed_rand_float(), Y_SCENE_SIZE*signed_rand_float(), rand_uniform(czmin, czmax));
				dirs[r] = (target_pt - pt).get_norm();
			}
			sort(dirs.begin(), dirs.end());

			for (unsigned r = 0; r < NRAYS; ++r) {
				if (dot_product(dirs[r], pt) >= 0.0) continue; // can get here when (-Z_SCENE_SIZE, Z_SCENE_SIZE) does not contain (czmin, czmax)
				point const end_pt(pt + dirs[r]*line_length);
				bool skip_this_line(0);

				// skip any lines that cross a global cube light because otherwise we would doubly add the contribution from that pos/dir
				for (vector<cube_light_source>::const_iterator i = global_cube_lights.begin(); i != global_cube_lights.end(); ++i) {
					if (check_line_clip(pt, end_pt, i->bounds.d)) {skip_this_line = 1; break;}
				}
				if (skip_this_line) continue;
				cast_light_ray(pt, end_pt, ray_wt, ray_wt, WHITE, line_length, -1, LIGHTING_SKY, 1);
				++start_rays;
			}
		}
	}
	for (vector<cube_light_source>::const_iterator i = global_cube_lights.begin(); i != global_cube_lights.end(); ++i) {
		if (data->num == 0) continue; // disabled
		if (data->verbose) cout << "Cube volume light source " << (i - global_cube_lights.begin()) << " of " << global_cube_lights.size() << endl;
		unsigned const num_rays(i->num_rays/data->num);
		point pt;

		for (unsigned p = 0; p < num_rays; ++p) {
			for (unsigned d = 0; d < 3; ++d) {
				pt[d] = rand_uniform(i->bounds.d[d][0], i->bounds.d[d][1]);
			}
			vector3d dir(signed_rand_vector_norm());
			dir.z = -fabs(dir.z); // make sure z is negative since this is supposed to be light from the sky
			point const end_pt(pt + dir*line_length);
			cast_light_ray(pt, end_pt, i->intensity, i->intensity, i->color, line_length, -1, LIGHTING_SKY, 1);
			++cube_start_rays;
		}
	}
	if (data->verbose) {
		cout << endl << "start rays: " << start_rays << "cube start rays: " << cube_start_rays << ", total rays: " << tot_rays
			 << ", hits: " << num_hits << ", cells touched: " << cells_touched << endl;
	}
	return 0;
}


void ray_trace_local_light_source(light_source const &ls, float line_length, unsigned num_rays) {

	point const &lpos(ls.get_center());
	colorRGBA lcolor(ls.get_color());
	float const ray_wt(1000.0*lcolor.alpha*ls.get_radius()/LOCAL_RAYS), r_inner(ls.get_r_inner());
	if (ray_wt == 0.0) return; // error?
	int init_cobj(-1);
	check_coll_line(lpos, lpos, init_cobj, -1, 1, 2); // find most opaque (max alpha) containing object
	assert(init_cobj < (int)coll_objects.size());
	
	for (unsigned n = 0; n < num_rays; ++n) {
		vector3d const dir(signed_rand_vector_spherical(1.0, 0).get_norm());
		float const weight(ray_wt*ls.get_dir_intensity(dir*-1));
		if (weight == 0.0) continue;
		point start_pt;

		if (init_cobj >= 0 && r_inner == 0.0) { // use a volume light source
			coll_obj cobj(coll_objects[init_cobj]);
			
			while (1) {
				for (unsigned d = 0; d < 3; ++d) { // generate a radom point within cobj's bounding cube
					start_pt[d] = rand_uniform(cobj.d[d][0], cobj.d[d][1]);
				}
				if (cobj.line_intersect(start_pt, start_pt)) break; // intersection (point contained)
			}
		}
		else {
			// move r_inner away from the light source
			// necessary for sources contained in more than one cobj (like the lamps in mapx)
			start_pt = lpos + dir*r_inner;
		}
		point const end_pt(start_pt + dir*line_length);
		cast_light_ray(start_pt, end_pt, weight, weight, lcolor, line_length, init_cobj, LIGHTING_LOCAL, 1);
	}
}


void *trace_ray_block_local(void *ptr) {

	assert(ptr);
	rt_data *data((rt_data *)ptr);
	cout << "Starting on thread " << data->ix << endl;
	assert(data->num > 0);
	if (data->is_thread) srand(data->rseed);
	float const scene_radius(get_scene_radius()), line_length(2.0*scene_radius);
	unsigned const num_rays(max(1U, LOCAL_RAYS/data->num));
	
	if (data->verbose) {
		cout << "Local light sources: " << light_sources.size() << ", progress: 0";
		cout.flush();
	}
	for (unsigned i = 0; i < light_sources.size(); ++i) {
		if (data->verbose) increment_printed_number(i);
		ray_trace_local_light_source(light_sources[i], line_length, num_rays);
	}
	if (data->verbose) cout << endl;
	return 0;
}


typedef void *(*ray_trace_func)(void *);
ray_trace_func const rt_funcs[NUM_LIGHTING_TYPES] = {trace_ray_block_sky, trace_ray_block_global, trace_ray_block_local};


void compute_ray_trace_lighting(unsigned ltype) {

	assert(ltype < NUM_LIGHTING_TYPES);

	if (read_light_files[ltype]) {
		lmap_manager.read_data_from_file(lighting_file[ltype], ltype);
	}
	else {
		if (ltype != LIGHTING_LOCAL) cout << X_SCENE_SIZE << " " << Y_SCENE_SIZE << " " << Z_SCENE_SIZE << " " << czmin << " " << czmax << endl;
		launch_threaded_job(NUM_THREADS, rt_funcs[ltype], 1);
	}
	if (write_light_files[ltype]) {
		lmap_manager.write_data_to_file(lighting_file[ltype], ltype);
	}
	lmap_manager.apply_light_scale(light_int_scale[ltype], ltype);
}


// lmap_manager_t


bool lmap_manager_t::read_data_from_file(char const *const fn, int ltype) {

	FILE *fp;
	assert(fn != NULL);
	if (!open_file(fp, fn, "lighting input", "rb")) return 0;
	cout << "Reading lighting file from " << fn << endl;
	unsigned data_size(0);
	size_t const sz_read(fread(&data_size, sizeof(unsigned), 1, fp));
	assert(sz_read == 1);

	if (data_size != vldata_alloc.size()) {
		cout << "Error: Lighting file " << fn << " data size of " << data_size
			 << " does not equal the expected size of " << vldata_alloc.size() << ". Ignoring file." << endl;
	}
	else {
		unsigned const sz(lmcell::get_dsz(ltype));

		for (vector<lmcell>::iterator i = vldata_alloc.begin(); i != vldata_alloc.end(); ++i) {
			size_t const nread(fread(i->get_offset(ltype), sizeof(float), sz, fp));
			assert(nread == sz);
		}
	}
	fclose(fp);
	return 1;
}


bool lmap_manager_t::write_data_to_file(char const *const fn, int ltype) const {

	FILE *fp;
	assert(fn != NULL);
	if (!open_file(fp, fn, "lighting output", "wb")) return 0;
	cout << "Writing lighting file to " << fn << endl;
	unsigned const data_size(vldata_alloc.size());
	size_t const sz_write(fwrite(&data_size, sizeof(unsigned), 1, fp));
	assert(sz_write == 1);
	unsigned const sz(lmcell::get_dsz(ltype));

	for (vector<lmcell>::const_iterator i = vldata_alloc.begin(); i != vldata_alloc.end(); ++i) { // const_iterator?
		size_t const nwrite(fwrite(i->get_offset(ltype), sizeof(float), sz, fp));
		assert(nwrite == sz);
	}
	fclose(fp);
	return 1;
}


void lmap_manager_t::apply_light_scale(float scale, int ltype) {

	assert(ltype < NUM_LIGHTING_TYPES);

	// apply global light scaling and normalize colors
	for (vector<lmcell>::iterator i = vldata_alloc.begin(); i != vldata_alloc.end(); ++i) {
		float *color(i->get_offset(ltype));

		if (ltype == LIGHTING_LOCAL) {
			UNROLL_3X(color[i_] = min(1.0f, (color[i_]*scale)););
			if (indir_light_exp != 0.0) {UNROLL_3X(color[i_] = pow(color[i_], indir_light_exp););} // gamma correction
		}
		else {
			if (color[3] == 0.0) continue;
			color[3] = CLIP_TO_01(color[3]*scale);
			if (indir_light_exp != 0.0) {color[3] = pow(color[3], indir_light_exp);} // gamma correction
			float const max_color(max(color[0], max(color[1], color[2])));
			if (max_color > 0.0) {UNROLL_3X(color[i_] /= max_color;)} // normalize color
		}
	}
}

