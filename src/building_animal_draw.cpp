// 3D World - Building Animal Drawing
// by Frank Gennari 3/5/23

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "shaders.h"
#include "city_model.h"

extern bool enable_depth_clamp, player_in_tunnel, player_in_mall;
extern int display_mode, animate2;
extern cube_t smap_light_clip_cube;
extern object_model_loader_t building_obj_model_loader;


class spider_draw_t {
	// could make mats[1] a subset of mats[0], but index logic is more complex
	rgeom_mat_t mats[2], web_mat; // mats: {high detail/67K verts, low detail shadow pass/4.2K verts}
	bool is_setup=0, had_inval_spider_warn=0;

	void add_eye(rgeom_mat_t &mat, point const &pos, float radius) {
		cube_t eye_bc(pos);
		eye_bc.expand_by(radius);
		mat.add_sphere_to_verts(eye_bc, DK_RED, 1); // low_detail=1
	}
	void init() {
		// generate spider geometry; centered at (0,0,0) with radius=1.0; head is in +X
		colorRGBA const color(BLACK);
		float const body_zval(-0.3), leg_radius(0.03);
		vector3d const sphere_radius(leg_radius, leg_radius, leg_radius);
		point const abdomen_center(-0.8, 0.0, body_zval);
		cube_t abdomen(abdomen_center), body(point(0.0, 0.0, body_zval));
		abdomen.expand_by(vector3d(0.50, 0.35, 0.35));
		body   .expand_by(vector3d(0.45, 0.30, 0.20));

		for (unsigned n = 0; n < 2; ++n) { // {high detail, low detail shadow pass}
			bool const low_detail(n == 1);
			unsigned const ndiv(get_rgeom_sphere_ndiv(low_detail)/2);
			unsigned cur_vert_pos(0);
			rgeom_mat_t &mat(mats[n]);
			mat.add_sphere_to_verts(abdomen, color, low_detail);
			mat.add_sphere_to_verts(body,    color, low_detail);

			if (!low_detail) { // red markings on back
				cube_t marking(abdomen_center + vector3d(-0.02, 0.0, 0.01));
				marking.expand_by(vector3d(0.48, 0.3, 0.35));
				mat.add_sphere_to_verts(marking, colorRGBA(0.5, 0.05, 0.0));
			}
			assign_tc_range(mat, cur_vert_pos, 0.0, 0.0, 0.0); // head and body aren't animated

			for (unsigned d = 0; d < 2; ++d) { // {left, right}
				float const d_sign(d ? -1.0 : 1.0);

				if (!low_detail) { // eyes and fangs are high detail only
					float const fang_radius(0.05);
					point const fang_top(0.44, 0.05*d_sign, body_zval-0.04), fang_bot(fang_top - vector3d(0.0, 0.0, 0.2));
					add_eye(mat, point(0.30, 0.080*d_sign, body_zval+0.14), 0.026);
					add_eye(mat, point(0.40, 0.045*d_sign, body_zval+0.08), 0.028);
					add_eye(mat, point(0.44, 0.020*d_sign, body_zval+0.04), 0.016);
					add_eye(mat, point(0.43, 0.055*d_sign, body_zval+0.03), 0.015);
					mat.add_sphere_to_verts(fang_top, vector3d(fang_radius, fang_radius, fang_radius), color, 1); // top of fang; low_detail=1
					mat.add_cylin_to_verts(fang_bot, fang_top, 0.0, fang_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv); // fang
					assign_tc_range(mat, cur_vert_pos, 0.0, 0.0, 0.0); // not animated
				}
				// add 4 pairs of legs
				for (unsigned l = 0; l < 4; ++l) {
					float const ts(l/4.0);
					point const joint(0.12*(l - 1.5), 0.26*d_sign, body_zval);
					point const knee (2.0*joint.x, 2.0*joint.y,  0.5);
					point const ankle(2.8*knee .x, 2.8*knee .y,  0.0);
					point const foot (3.5*knee .x, 3.5*knee .y, -1.0);
					float const joint_tt(0.0*d_sign), knee_tt(0.3*d_sign), ankle_tt(0.7*d_sign), foot_tt(1.0*d_sign);
					if (!low_detail) {mat.add_sphere_to_verts(joint, sphere_radius, color, 1);} // round body joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, joint_tt, joint_tt);
					mat.add_cylin_to_verts(joint, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, joint_tt, knee_tt);
					if (!low_detail) {mat.add_sphere_to_verts(knee, sphere_radius, color, 1);} // round knee joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, knee_tt, knee_tt);
					mat.add_cylin_to_verts(ankle, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, ankle_tt, knee_tt);
					if (!low_detail) {mat.add_sphere_to_verts(ankle, sphere_radius, color, 1);} // round ankle joint; high detail only; low_detail=1
					assign_tc_range(mat, cur_vert_pos, ts, ankle_tt, ankle_tt);
					mat.add_cylin_to_verts(foot,  ankle, 0.1*leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
					assign_tc_range(mat, cur_vert_pos, ts, foot_tt, ankle_tt);
				} // for l
			} // for d
			mat.create_vbo_inner();
			mat.clear_vectors(1); // free_memory=1: vector data no longer needed
		} // for n
		is_setup = 1;
	}
	void assign_tc_range(rgeom_mat_t &mat, unsigned &cur_vert_pos, float ts, float tt_lo, float tt_hi) {
		// tcs are used for animation:
		// ts is the position of the leg from front to back: {0.0, 0.25, 0.5, 0.75}
		// tt is the joint index: negative for left, positive for right; 0.0 for body joint, 0.5 for knee, 1.0 for foot
		for (auto i = mat.itri_verts.begin()+cur_vert_pos; i != mat.itri_verts.end(); ++i) {
			i->t[0] = ts;
			i->t[1] = i->t[1]*tt_hi + (1.0 - i->t[1])*tt_lo; // existing t[1] is 0.0 for low end and 1.0 for high end
		}
		cur_vert_pos = mat.itri_verts.size();
	}
public:
	void clear() {mats[0].clear(); mats[1].clear(); is_setup = 0;}

	void draw(vect_spider_t const &spiders, shader_t &s, building_t const &building, occlusion_checker_noncity_t const &oc,
		vector3d const &xlate, bool shadow_only, bool reflection_pass, bool check_clip_cube)
	{
		if (spiders.empty()) return; // nothing to draw
		int anim_time_loc(-1);
		point const camera_bs(camera_pdu.pos - xlate);
		bool const check_occlusion(display_mode & 0x08), low_detail(shadow_only || reflection_pass);
		bool any_drawn(0);
		rgeom_mat_t &mat(mats[low_detail]);

		for (spider_t const &S : spiders) { // future work: use instancing
			if (!shadow_only && S.on_web) { // draw spider webs
				cube_t strand(S.pos - 0.3*S.radius*S.upv); // end of abdomen (matches body_zval)
				set_cube_zvals(strand, (S.pos.z + 1.3*S.radius), S.web_start_zval);
				strand.expand_by_xy(0.02*S.radius);

				if (strand.z1() < strand.z2() && camera_pdu.cube_visible(strand + xlate)) {
					web_mat.add_vcylin_to_verts(strand, WHITE, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 16); // ndiv=16
				}
			}
			if (!player_in_mall && S.in_tank)        continue; // player not in mall, likely not visible
			if (shadow_only     && S.in_tank)        continue; // too small to cast a shadow
			if (shadow_only && S.shadow_non_visible) continue; // shadow not visible to the camera (player)
			cube_t const bcube(S.get_bcube());
			if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
			if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
			if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass, 0, 0, 1)) continue; // inc_extra_occluders=1
			
			if (S.dir == zero_vector || S.upv == zero_vector) {
				if (!had_inval_spider_warn) {cout << "Error: Invalid spider: " << TXTS(S.pos) << TXTS(S.dir) << TXTS(S.upv) << endl;} // print only once
				had_inval_spider_warn = 1;
				continue; // seems like this can occasionally happen with too many objects and low framerate; maybe FP error; make it nonfatal
			}
			if (!any_drawn) { // setup shaders
				if (!is_setup) {init();}
				mat.vao_setup(shadow_only);
				s.set_specular(0.5, 80.0);
				select_no_texture();
				int const animation_id = ANIM_ID_SPIDER; // custom spider animation; not using animation_state_t here
				s.add_uniform_int("animation_id", animation_id);
				s.add_uniform_float("animation_scale",    1.0); // not using a model, nominal size is 1.0
				s.add_uniform_float("model_delta_height", 1.0); // not using a model, nominal size is 1.0
				s.add_uniform_float("bump_map_mag",       0.0);
				mat.pre_draw(shadow_only);
				anim_time_loc = s.get_uniform_loc("animation_time");
				any_drawn = 1;
			}
			s.set_uniform_float(anim_time_loc, S.anim_time);
			fgPushMatrix();
			translate_to(S.pos);
			vector3d const dir(S.dir.get_norm()), upv(S.upv.get_norm()); // normalize, just in case
			vector3d const right(cross_product(S.upv, S.dir).get_norm());
			assert(right != zero_vector);
			xform_matrix xm;
			float *m(xm.get_ptr());
			m[0] = dir.x; // x
			m[1] = dir.y;
			m[2] = dir.z;
			m[4] = right.x; // y
			m[5] = right.y;
			m[6] = right.z;
			m[8] = upv.x; // z
			m[9] = upv.y;
			m[10]= upv.z;
			fgMultMatrix(xm);
			uniform_scale(S.radius);
			if (S.squished) {fgScale(1.0, 1.0, 0.2);} // flatten it if squished
			check_mvm_update();
			// Note: we could use hardware instancing here if the MVM is used as a uniform in the shader, but we also need to set per-instance anim_time
			mat.draw_geom();
			fgPopMatrix();
		} // for S
		if (any_drawn) { // reset state
			check_mvm_update(); // make sure to reset MVM
			s.add_uniform_float("bump_map_mag",   1.0);
			s.add_uniform_float("animation_time", 0.0); // reset animation time
			s.add_uniform_int  ("animation_id",   0); // clear animation
			s.clear_specular();
			indexed_vao_manager_with_shadow_t::post_render();
		}
		if (!web_mat.empty()) {
			select_no_texture();
			tid_nm_pair_dstate_t state(s);
			web_mat.upload_draw_and_clear(state);
		}
	}
};
spider_draw_t spider_draw;

class insect_draw_t {
	rgeom_mat_t fly_mat, roach_mat;
	float timebase=0.0;
	vector<vert_norm_comp> wing_verts; // temp used in draw calls for flies
	vector<unsigned> to_draw[NUM_INSECT_TYPES]; // temp state used when drawing

	void init_fly() { // generate fly geometry, just a simple sphere for now
		if (fly_mat.num_verts > 0) return; // already setup
		bool const low_detail = 1;
		colorRGBA const color(BLACK);
		float const body_zval(0.0), leg_radius(0.02);
		cube_t thorax(point(0.2, 0.0, body_zval)), abdomen(point(-0.3, 0.0, body_zval-0.1)), head(point(0.8, 0.0, body_zval));
		thorax .expand_by(vector3d(0.50, 0.33, 0.28));
		abdomen.expand_by(vector3d(0.70, 0.22, 0.22));
		head   .expand_by(vector3d(0.18, 0.26, 0.24));
		fly_mat.add_sphere_to_verts(thorax,  color, low_detail);
		fly_mat.add_sphere_to_verts(abdomen, color, low_detail);
		fly_mat.add_sphere_to_verts(head,    color, low_detail);

		for (unsigned d = 0; d < 2; ++d) { // {left, right}
			float const d_sign(d ? -1.0 : 1.0);
			point const eye_pos(0.90, 0.12*d_sign, 0.05);
			cube_t eye(eye_pos);
			eye.expand_by(vector3d(0.08, 0.14, 0.14));
			fly_mat.add_sphere_to_verts(eye, colorRGBA(0.5, 0.1, 0.0), low_detail); // dark red-orange
			// add 3 pairs of legs; don't need to draw the leg joints because flies are so small
			unsigned const ndiv = 8;

			for (unsigned l = 0; l < 3; ++l) {
				point const joint(0.12*(l - 1.5), 0.26*d_sign, body_zval);
				point const knee (2.0*joint.x, 2.0*joint.y,  0.25);
				point const ankle(1.5*knee .x, 1.5*knee .y, -0.20);
				point const foot (3.0*knee .x, 3.0*knee .y, -0.70);
				fly_mat.add_cylin_to_verts(joint, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
				fly_mat.add_cylin_to_verts(ankle, knee, leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
				fly_mat.add_cylin_to_verts(foot,  ankle, 0.1*leg_radius, leg_radius, color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
			} // for l
		} // for d
		fly_mat.create_vbo_inner();
		fly_mat.clear_vectors(1); // free_memory=1: vector data no longer needed
	}
	void init_roach() {
		if (roach_mat.num_verts > 0) return; // already setup
		bool const low_detail = 0;
		colorRGBA const color(0.2, 0.1, 0.05, 1.0); // darker brown
		cube_t body; // centered on (0,0,0)
		body.expand_by(vector3d(1.2, 0.5, 0.2));
		roach_mat.add_sphere_to_verts(body, color, low_detail, -plus_z); // skip bottom
		// cockroaches are fast and run when disturbed, and their legs are under them, so they don't really need to be drawn/animated
		roach_mat.create_vbo_inner();
		roach_mat.clear_vectors(1); // free_memory=1: vector data no longer needed
	}
	void draw_insect_list(vect_insect_t const &insects, shader_t &s, rgeom_mat_t &mat, vector3d const &xlate,
		vector<unsigned> const &ixs, unsigned obj_model_ix=NUM_OBJ_MODELS) const
	{
		if (ixs.empty()) return;
		bool const use_model(obj_model_ix < NUM_OBJ_MODELS && building_obj_model_loader.is_model_valid(obj_model_ix));

		if (!use_model) { // no model, setup material for drawing
			mat.vao_setup(0); // shadow_only=0
			mat.pre_draw (0); // shadow_only=0
		}
		for (unsigned ix : ixs) { // for ix
			assert(ix < insects.size());
			insect_t const &i(insects[ix]);
			
			if (use_model) {
				cube_t const bcube(i.get_bcube_with_dir());
				building_obj_model_loader.draw_model(s, cube_bot_center(bcube), bcube, i.dir, LT_BROWN, xlate, obj_model_ix, 0); // shadow_only=0
			}
			else { // draw as untextured sphere geometry
				//if (i.has_target) {s.set_color_e(i.target_player ? RED : GREEN);} // debug visualization for flies
				fgPushMatrix();
				translate_to(i.pos);
				rotate_from_v2v(i.get_orient(), plus_x); // rotate around Z axis
				uniform_scale(i.radius);
				check_mvm_update();
				mat.draw_geom(); // use hardware instancing?
				fgPopMatrix();
				//if (i.has_target) {s.clear_color_e();}
			}
		} // for ix
		check_mvm_update();
	}
public:
	void clear() {fly_mat.clear(); roach_mat.clear();}

	void draw(vect_insect_t const &insects, shader_t &s, building_t const &building, occlusion_checker_noncity_t const &oc, vector3d const &xlate, bool reflection_pass) {
		if (insects.empty()) return; // nothing to draw
		point const camera_bs(camera_pdu.pos - xlate);
		bool const check_occlusion(display_mode & 0x08);//, low_detail(shadow_only || reflection_pass);
		float const draw_dist_scale = 500.0;
		bool any_drawn(0);
		wing_verts.clear();
		if (animate2) {timebase = tfticks;}

		for (auto i = insects.begin(); i != insects.end(); ++i) { // future work: use instancing
			if (!dist_less_than(i->pos, camera_bs, draw_dist_scale*i->radius)) continue; // too far
			cube_t const bcube(i->get_bcube());
			if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
			if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass)) continue;
			assert(i->type < NUM_INSECT_TYPES);
			to_draw[i->type].push_back(i - insects.begin());
			any_drawn = 1;
		}
		if (!any_drawn) return;
		select_no_texture();
		s.add_uniform_float("bump_map_mag", 0.0); // no normal maps

		if (!to_draw[INSECT_TYPE_FLY].empty()) { // draw flies
			init_fly();
			s.set_specular(0.5, 80.0);
			if (!enable_depth_clamp) {glEnable(GL_DEPTH_CLAMP);} // make sure depth clamp is enabled so that insects are drawn when very close
			draw_insect_list(insects, s, fly_mat, xlate, to_draw[INSECT_TYPE_FLY]);

			for (unsigned ix : to_draw[INSECT_TYPE_FLY]) { // draw the wings
				insect_t const &i(insects[ix]);
				if (!dist_less_than(i.pos, camera_bs, 0.25*draw_dist_scale*i.radius)) continue; // too far to draw wings
				vector3d const orient(i.get_orient());
				vector3d const side_dir(cross_product(orient, plus_z)); // should be normalized
				norm_comp const normal(plus_z); // use actual normal?
				float const lift_amt(0.5 + 0.5*sin(4.0*(i.anim_time + timebase))); // add in global time so that wings still flap when hovering

				for (unsigned d = 0; d < 2; ++d) { // {left, right}
					float const d_sign(d ? -1.0 : 1.0);
					point v[3]; // wing triangle verts, in local coordinate space of fly model
					v[0] =  0.25*orient + 0.10*d_sign*side_dir + 0.3*plus_z; // back connect point
					v[1] = -0.30*orient - 0.02*d_sign*side_dir + 0.3*plus_z; // near center of body
					v[2] = -1.50*orient + (1.0 - 0.7*lift_amt)*d_sign*side_dir + (0.3 + 0.7*lift_amt)*plus_z; // tip
					UNROLL_3X(wing_verts.emplace_back((i.pos + i.radius*v[i_]), normal););
				} // for d
			} // for ix
			to_draw[INSECT_TYPE_FLY].clear();
			indexed_vao_manager_with_shadow_t::post_render(); // unbind VBO/VAO

			if (!wing_verts.empty()) {
				glDisable(GL_CULL_FACE); // wings are two sided
				enable_blend();
				s.set_cur_color(colorRGBA(1.0, 1.0, 1.0, 0.25)); // transparent white
				draw_verts(wing_verts, GL_TRIANGLES);
				s.set_cur_color(WHITE);
				disable_blend();
				glEnable(GL_CULL_FACE);
			}
			if (!enable_depth_clamp) {glDisable(GL_DEPTH_CLAMP);}
			s.clear_specular();
		}
		if (!to_draw[INSECT_TYPE_ROACH].empty()) { // draw cockroaches
			bool const draw_as_model(building_obj_model_loader.is_model_valid(INSECT_TYPE_ROACH));
			if (!draw_as_model) {init_roach();} // only setup if not drawing the 3D model
			if (!draw_as_model) {s.set_specular(0.35, 40.0);}
			draw_insect_list(insects, s, roach_mat, xlate, to_draw[INSECT_TYPE_ROACH], OBJ_MODEL_ROACH);
			if (!draw_as_model) {s.clear_specular();}
			to_draw[INSECT_TYPE_ROACH].clear();
			indexed_vao_manager_with_shadow_t::post_render();
		}
		// reset state
		check_mvm_update(); // make sure to reset MVM
		s.add_uniform_float("bump_map_mag", 1.0);
	}
};
insect_draw_t insect_draw;

// Note: similar to the functions in Tree.cpp, but pushes back rather than assigning, and step is hard-coded to 1
void add_cylin_indices_tris(vector<unsigned> &idata, unsigned ndiv, unsigned ix_start) {
	for (unsigned S = 0; S < ndiv; ++S) {
		bool const last_edge(S == ndiv-1);
		unsigned const ix0(ix_start + S), ixs[4] = {0, ndiv, (last_edge ? 1 : 1+ndiv), (last_edge ? 1-ndiv : 1)};
		for (unsigned i = 0; i < 6; ++i) {idata.push_back(ix0 + ixs[quad_to_tris_ixs[i]]);}
	}
}
// used for snakes and vases
void draw_segment(rgeom_mat_t &mat, point const &p1, point const &p2, float radius1, float radius2,
	float seg_ix, float tscale_x, float tscale_y, color_wrapper const &cw, unsigned ndiv, unsigned &data_pos)
{
	point const ce[2] = {p1, p2};
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius1, radius2, ndiv, v12, NULL, 0.0, 1.0, 2)); // force_dim=2
	float const ndiv_inv(1.0/ndiv);
	bool const is_first(seg_ix == 0.0);

	for (unsigned j = !is_first; j < 2; ++j) {
		float const ty(tscale_y*(seg_ix + j));

		for (unsigned S = 0; S < ndiv; ++S) {
			float const tx(tscale_x*fabs(S*ndiv_inv - 0.5f));
			vector3d const n(0.5f*(vpn.n[S] + vpn.n[(S+ndiv-1)%ndiv])); // average face normals to get vert normals, don't need to normalize
			mat.itri_verts.emplace_back(vpn.p[(S<<1)+j], n, tx, ty, cw);
		}
	}
	add_cylin_indices_tris(mat.indices, ndiv, data_pos); // create index data
	data_pos += ndiv;
}

class snake_draw_t {
	rgeom_mat_t skin_mats[2], untex_mat;

	void draw_snake(snake_t const &S, bool shadow_only, bool reflection_pass, bool is_distant) {
		bool const low_detail(shadow_only || reflection_pass || is_distant);
		unsigned const ndiv(get_rgeom_sphere_ndiv(low_detail)/2);
		float const tscale = 1.0; // must tune this once our snake is textured
		colorRGBA color(S.color);
		//if (S.stuck_counter > 0.0) {color = blend_color(RED, S.color, min(1.0, 0.0025*S.stuck_counter), 0);} // debugging: turn red over time when stuck
		// draw head
		float const head_hheight(0.6*S.radius), head_hlen(1.6*S.radius), head_hwidth(S.radius);
		vector3d const &dir(S.last_valid_dir);
		vector3d const head_size(head_hlen, head_hwidth, head_hheight); // length in X, max radius in Y, flattened in Z
		point const head_pos(S.get_head_pos()), head_center(head_pos + vector3d(0,0,head_hheight));
		rgeom_mat_t &skin_mat(skin_mats[S.id & 1]); // use the ID's LSB to select which of the two skin textures will be used
		unsigned const head_verts_start(skin_mat.itri_verts.size());
		skin_mat.add_sphere_to_verts(head_center, head_size, color, low_detail);
		float const rot_angle(-atan2(dir.y, dir.x)); // rotate dir into +X
		rotate_verts(skin_mat.itri_verts, plus_z, rot_angle, head_pos, head_verts_start);
		// draw segments
		float const zscale = 0.85;
		color_wrapper const cw(color);
		unsigned const body_verts_start(skin_mat.itri_verts.size()), rattle_verts_start(untex_mat.itri_verts.size());
		unsigned data_pos(body_verts_start);
		vector3d prev_v12;

		for (auto s = S.segments.begin(); s != S.segments.end(); ++s) {
			unsigned const seg_ix(s - S.segments.begin());
			bool const is_first(seg_ix == 0), is_tail(s+1 == S.segments.end());
			float const radius1(S.get_seg_radius(seg_ix)), radius2(S.get_seg_radius(seg_ix+1));
			point const seg_start(is_first ? *s : 0.5*(*s + *(s-1))); // midpoint between this segment and the last
			point const seg_end  (is_tail  ? (*s + 2.0*(*s - seg_start)) : 0.5*(*s + *(s+1))); // midpoint between this segment and the next; tail extends further back
			bool const add_rattle(S.has_rattle && s+1 == S.segments.end() && !shadow_only);
			point const ce[2] = {(seg_start + vector3d(0,0,radius1)), (seg_end + vector3d(0,0,radius2))};

			if (add_rattle) { // add a rattle on the last tail segment if not the shadow pass
				unsigned const num_segs(5 + (S.id % 4)); // 5-8
				vector3d const seg_delta(ce[1] - ce[0]);
				colorRGBA const rattle_color(1.0, 0.8, 0.6); // yellow-ish
				vector3d const step(seg_delta/num_segs);
				float const rattle_radius(0.9*step.mag()); // add some overlap between segments
				
				for (unsigned i = 0; i < num_segs; ++i) {
					point rattle_pos(ce[0] + (i+0.5)*step);
					float const seg_radius(rattle_radius*(1.0 - 0.5*i/float(num_segs)));
					rattle_pos.z = head_pos.z + seg_radius; // place exactly on the floor
					untex_mat.add_sphere_to_verts(rattle_pos, seg_radius, rattle_color, 1); // low_detail=1
				}
			}
			else { // skip tail if adding a rattle
				vector3d const delta(ce[1] - ce[0]), v12(delta.get_norm());

				if (!is_first && !is_tail && dot_product(v12, prev_v12) < 0.75) { // sharp bend, draw as two cylinders
					vector3d const v12_avg(0.5f*(v12 + prev_v12)); // use the average vector for a more gradual transition
					float const r_mid(0.5*(radius1 + radius2));
					point const seg_center(ce[0] + v12_avg*(0.5*delta.mag()));
					draw_segment(skin_mat, ce[0], seg_center, radius1, r_mid, seg_ix+0.0, 2.0, tscale, cw, ndiv, data_pos);
					draw_segment(skin_mat, seg_center, ce[1], r_mid, radius2, seg_ix+0.5, 2.0, tscale, cw, ndiv, data_pos);
				}
				else {
					draw_segment(skin_mat, ce[0], ce[1], radius1, radius2, seg_ix, 2.0, tscale, cw, ndiv, data_pos);
				}
				prev_v12 = v12;
			}
		} // for s
		if (zscale != 1.0) { // flatten a bit in Z
			for (auto i = skin_mat .itri_verts.begin()+body_verts_start  ; i != skin_mat.itri_verts .end(); ++i) {i->v.z =      zscale*(i->v.z - head_pos.z) + head_pos.z;}
			for (auto i = untex_mat.itri_verts.begin()+rattle_verts_start; i != untex_mat.itri_verts.end(); ++i) {i->v.z = 0.7f*zscale*(i->v.z - head_pos.z) + head_pos.z;}
		}
		if (low_detail) return; // no eyes or tongue in low detail mode
		// add eyes to head
		vector3d const side_dir(cross_product(dir, plus_z));
		float const eye_extend_fwd(0.9*head_hlen/SQRT2), eye_extend_side(0.86*head_hwidth/SQRT2);

		for (unsigned d = 0; d < 2; ++d) {
			point eye_pos(head_center);
			eye_pos += eye_extend_fwd*dir; // move forward
			eye_pos += ((d ? -1.0 : 1.0)*eye_extend_side)*side_dir; // move to the side
			untex_mat.add_sphere_to_verts(eye_pos, 0.25*S.radius, BLACK, 1); // low_detail=1
		}
		if (fract(0.015f*S.anim_time) > 0.8) { // occasionally stick tongue out
			colorRGBA const tongue_color((S.id & 2) ? BLACK : colorRGBA(1.0, 0.6, 0.4)); // black or pink
			color_wrapper const tip_color(BLACK);
			float const radius(0.04*head_hlen), fwd_dist(0.9*head_hlen), split_len(0.6*head_hlen), tip_len(0.5*head_hlen);
			
			for (unsigned d = 0; d < 2; ++d) { // draw each half as a cylinder + cone
				point const start(head_center + fwd_dist*dir + (d ? -1.0 : 1.0)*0.75*radius*side_dir);
				point const split(start + split_len*dir);
				point const tip(split + tip_len*(dir + (d ? -1.0 : 1.0)*0.67*side_dir));
				untex_mat.add_cylin_to_verts(start, split, radius, radius, tongue_color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);
				unsigned const tip_verts_start(untex_mat.itri_verts.size());
				untex_mat.add_cylin_to_verts(split, tip,   radius,    0.0, tongue_color, 0, 0, 0, 0, 1.0, 1.0, 0, ndiv);

				if (tongue_color != BLACK) { // recolor the tip, every second vertex
					for (unsigned i = tip_verts_start+1; i < untex_mat.itri_verts.size(); i += 2) {untex_mat.itri_verts[i].copy_color(tip_color);}
				}
				untex_mat.add_sphere_to_verts(split, radius, tongue_color, 1);
			} // for d
		}
	}
public:
	void draw(vect_snake_t const &snakes, vect_snake_t const &pet_snakes, shader_t &s, building_t const &building, occlusion_checker_noncity_t const &oc,
		vector3d const &xlate, bool shadow_only, bool reflection_pass, bool check_clip_cube)
	{
		if (snakes.empty() && pet_snakes.empty()) return; // nothing to draw
		//highres_timer_t timer("Draw Snakes");
		point const camera_bs(camera_pdu.pos - xlate);
		bool const check_occlusion(display_mode & 0x08), draw_pet_snakes(!shadow_only && player_in_mall); // skip pet snakes in the shadow pass
		bool any_drawn(0);

		for (unsigned d = 0; d < (draw_pet_snakes ? 2U : 1U); ++d) {
			for (snake_t const &S : (d ? pet_snakes : snakes)) {
				if (shadow_only && S.shadow_non_visible) continue; // shadow not visible to the camera (player)
				cube_t const bcube(S.get_bcube());
				if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
				if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
				if (check_occlusion && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass, 0, 0, (d == 1))) continue; // inc_extra_occluders=pet snakes
				bool const is_distant(!dist_less_than(camera_bs, S.pos, 6.0*S.length));
				draw_snake(S, shadow_only, reflection_pass, is_distant);
				any_drawn = 1;
			} // for S
		} // for d
		if (!any_drawn) return;
		tid_nm_pair_dstate_t state(s);
		s.add_uniform_float("bump_map_mag", 0.0);
		// draw the skin material
		s.set_specular(0.25, 50.0);

		for (unsigned d = 0; d < 2; ++d) {
			rgeom_mat_t &skin_mat(skin_mats[d]);
			if (skin_mat.tex.tid < 0) {skin_mat.tex = tid_nm_pair_t(get_texture_by_name(d ? "interiors/snakeskin2.jpg" : "interiors/snakeskin.jpg"), 0.0, 1);}
			skin_mat.upload_draw_and_clear(state);
		}
		// draw the eyes and rattle
		s.set_specular(0.75, 80.0);
		untex_mat.upload_draw_and_clear(state);
		s.add_uniform_float("bump_map_mag", 1.0);
		s.clear_specular();
	}
};
snake_draw_t snake_draw;

void building_room_geom_t::draw_animals(shader_t &s, building_t const &building, occlusion_checker_noncity_t const &oc, vector3d const &xlate,
	point const &camera_bs, bool shadow_only, bool reflection_pass, bool check_clip_cube) const
{
	// would it help to pass in wall occluders for pet stores?
	if ((!rats.empty() || !sewer_rats.empty() || !pet_rats.empty()) && building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) {
		bool const enable_animations(!shadow_only); // can't see the animation in the shadow pass
		animation_state_t anim_state(enable_animations, ANIM_ID_RAT);
		bool rat_drawn(0);

		for (rat_t const &rat : rats) {
			if (shadow_only && rat.shadow_non_visible) continue; // shadow not visible to the camera (player)
			if (rat.in_drain_amt == 1.0) continue; // fully inside drain
			cube_t bcube(rat.get_bcube());
			if (rat.in_drain_amt > 0.0) {bcube.translate_dim(2, -rat.in_drain_amt*rat.height);} // sink into drain
			if (check_clip_cube && !smap_light_clip_cube.intersects(bcube + xlate)) continue; // shadow map clip cube test: fast and high rejection ratio, do this first
			if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
			if ((display_mode & 0x08) && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass)) continue;
			point const pos(bcube.get_cube_center());
			anim_state.anim_time = rat.anim_time;
			cube_t rat_bcube(rat.get_bcube_with_dir());
			// translate flooring thickness up so that feet don't clip through flooring and rugs as much, though now rats will slightly float above normal floors
			rat_bcube.translate_dim(2, building.get_flooring_thick());
			vector3d dir(rat.dir);

			if (rat.in_drain_amt > 0.0) { // partially inside drain
				rat_bcube.translate_dim(2, -rat.in_drain_amt*rat.height); // sink into drain
				bool const exiting(rat.fear == 0.0);
				dir.z += rat.in_drain_amt*(exiting ? 1.0 : -1.0);
				dir.normalize();
			}
			colorRGBA const color(rat_color); // make the rat's fur darker
			//colorRGBA const color(blend_color(RED, WHITE, rat.fear, 0)); // used for debugging fear
			//colorRGBA const color(blend_color(RED, WHITE, rat.attacking, 0));
			building_obj_model_loader.draw_model(s, pos, rat_bcube, dir, color, xlate, OBJ_MODEL_RAT, shadow_only, 0, &anim_state, 0, 0, 0, rat.dead); // upside down if dead

			if (rat.attacking) { // draw red glowing eyes
				s.set_color_e(colorRGBA(0.5, 0.0, 0.0, 1.0)); // light emissive red
				s.set_cur_color(RED);
				select_no_texture();
				anim_state.clear_animation_id(s); // clear animations
				point eyes_center(pos + vector3d(0.0, 0.0, 0.09*rat.height) + 0.85*rat.get_hlength()*dir);
				vector3d const eye_sep_dir(0.21*rat.hwidth*cross_product(dir, plus_z).get_norm());

				for (unsigned d = 0; d < 2; ++d) { // draw left and right eye, untextured
					draw_sphere_vbo((eyes_center + (d ? 1.0 : -1.0)*eye_sep_dir), 0.05*rat.height, 16, 0);
				}
				s.set_color_e(BLACK);
			}
			rat_drawn = 1;
		} // for rat
		if (!shadow_only) { // draw sewer and pet rats; not in the shadow pass
			vect_rat_t const *to_draw(nullptr);
			if      (player_in_tunnel) {to_draw = &sewer_rats;}
			else if (player_in_mall  ) {to_draw = &pet_rats  ;}

			if (to_draw) {
				colorRGBA const pet_rat_colors[4] = {WHITE, LT_BROWN, GRAY, BLACK};

				for (rat_t const &rat : *to_draw) {
					cube_t const bcube(rat.get_bcube());
					if (!camera_pdu.cube_visible(bcube + xlate)) continue; // VFC
					if ((display_mode & 0x08) && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass, 0, 0, 1)) continue; // inc_extra_occluders=1
					anim_state.anim_time = rat.anim_time;
					colorRGBA const color(player_in_mall ? pet_rat_colors[rat.id&3] : rat_color);
					building_obj_model_loader.draw_model(s, bcube.get_cube_center(), rat.get_bcube_with_dir(), rat.dir, color, xlate, OBJ_MODEL_RAT, shadow_only, 0, &anim_state);
					rat_drawn = 1;
				}
			}
		}
		anim_state.clear_animation_id(s); // clear animations
		bind_default_flat_normal_map();
		if (rat_drawn) {check_mvm_update();} // needed after popping model transform matrix
	} // end rats drawing
	if (!shadow_only && !pet_birds.empty() && building_obj_model_loader.is_model_valid(OBJ_MODEL_BIRD_ANIM)) {
		bool const enable_animations(1);
		animation_state_t anim_state(enable_animations, ANIM_ID_SKELETAL, 0.0, BIRD_STATE_STANDING);
		bool bird_drawn(0);

		for (pet_bird_t const &bird : pet_birds) {
			if (!camera_pdu.sphere_visible_test((bird.pos + xlate), bird.radius)) continue; // VFC
			cube_t bcube(bird.pos);
			bcube.expand_by(bird.radius);
			if ((display_mode & 0x08) && building.check_obj_occluded(bcube, camera_bs, oc, reflection_pass, 0, 0, 1)) continue; // inc_extra_occluders=1
			anim_state.anim_time = 0.02*bird.anim_time/SKELETAL_ANIM_TIME_CONST;
			building_obj_model_loader.draw_model(s, bird.pos, bcube, bird.dir, bird.color, xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);
			bird_drawn = 1;
		} // for bird
		anim_state.clear_animation_id(s); // clear animations
		bind_default_flat_normal_map();
		if (bird_drawn) {check_mvm_update();} // needed after popping model transform matrix
	} // end birds drawing
	spider_draw.draw(spiders,            s, building, oc, xlate, shadow_only, reflection_pass, check_clip_cube);
	snake_draw .draw(snakes, pet_snakes, s, building, oc, xlate, shadow_only, reflection_pass, check_clip_cube);
	if (!shadow_only) {insect_draw.draw(insects, s, building, oc, xlate, reflection_pass);} // insects are too small to cast shadows
	// draw sewer spiders
	if (player_in_tunnel && !shadow_only) {spider_draw.draw(sewer_spiders, s, building, oc, xlate, shadow_only, reflection_pass, check_clip_cube);}
}


