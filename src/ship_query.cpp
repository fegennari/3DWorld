// 3D World - Ship models for universe mode
// by Frank Gennari
// 8/25/05

#include "ship.h"
#include "ship_util.h"
#include "explosion.h"
#include "obj_sort.h"


bool const EXPLODE_LIGHTING = 1;


extern int display_mode;
extern float uobj_rmax, urm_ship, urm_static, urm_proj, urm_nstat;
extern vector<cached_obj> ships[], all_ships, stat_objs, coll_proj, decoys, c_uobjs;
extern vector<usw_ray> t_wrays;
extern vector<us_weapon> us_weapons;



// what about objects created this frame that aren't sorted?
unsigned binary_search_pos(vector<cached_obj> const &objs, point const &pos) { // returns the index before

	unsigned const size(objs.size());
	unsigned start, end;

	for (start = 0, end = size; ((end - start) > 1);) {
		unsigned const mid((start + end) >> 1);
		assert(start <= end && mid < size);
		if (pos.x > objs[mid].pos.x) {start = mid;} else {end = mid;}
	}
	return start;
}


void line_intersect_fo_vector(line_int_data &li_data, vector<cached_obj> const &objs, free_obj *&fobj, float urm, bool find_ships) {

	unsigned const nobjs(objs.size());
	if (nobjs == 0) return;
	float const line_radius(li_data.line_radius);
	urm += line_radius;
	bool const sign(li_data.dir.x > 0);
	int const ie(sign ? nobjs+1 : 0), di(sign ? 1 : -1);
	point start2(li_data.start);
	float const st_val(li_data.start.x), dmax(fabs(li_data.end.x - st_val) + 1.2*urm); // 2.0*urm?
	start2.x -= 1.01*di*urm;
	unsigned bad_flags(OBJ_FLAGS_BAD_); // Note: Bad (dying) objects can still get in the way
	if (!li_data.even_ncoll) bad_flags |= OBJ_FLAGS_NCOL;
	if (!find_ships)         bad_flags |= OBJ_FLAGS_SHIP;
	unsigned const six(binary_search_pos(objs, start2)); // could store the sort index in the object?
	vector<uobject const *> *sobjs(li_data.sobjs);
	vector3d const v_line(li_data.start, li_data.end);
	float t_val; // unused

	for (int i = six; i+1 != ie; i += di) {
		cached_obj const &obj(objs[i]);
		if (obj.flags & bad_flags) continue; // already destroyed or no collisions
		assert(obj.obj != NULL);
		if (obj.obj == li_data.curr || obj.obj == li_data.ignore_obj) continue; // don't hit yourself or ignore_obj
		point const &pos(obj.pos);

		// since we're using start2, not start, have to make sure we're comparing in the correct direction
		// also, objs created this frame aren't sorted, so can't break on them
		if (!(obj.flags & OBJ_FLAGS_NEW_) && ((st_val > pos.x) ^ sign)) { // move up?
			if (fabs(st_val - pos.x) > dmax) break; // critical performance improvement
		}
		float const radius(obj.radius + line_radius), rdist(radius + li_data.length), dist_sq(p2p_dist_sq(li_data.start, pos));
		if (dist_sq > rdist*rdist || (fobj != NULL && sobjs == NULL && dist_sq >= li_data.dist)) continue;

		// check_parent: 0 = disabled, 1 = projectiles only, 2 = projectiles + fighters
		if (li_data.check_parent && (li_data.check_parent == 2 || (obj.flags & OBJ_FLAGS_PROJ)) &&
			obj.obj->get_root_parent() == li_data.curr)
		{
			continue; // don't hit your own shot/fighter
		}
		if (!sphere_test_comp(li_data.start, pos, v_line, radius*radius, t_val))                 continue;
		if (li_data.visible_only && (obj.flags & OBJ_FLAGS_SHIP) && obj.obj->visibility() < 0.1) continue; // cache miss, rarely fails

		if (line_radius == 0.0 || !li_data.use_lpos) {
			if (!obj.obj->line_int_obj(li_data.start, li_data.end)) continue; // skip this check for thick lines
		}
		else { // thick lines, used for shadow calculations
			vector3d_d const test_dir((li_data.lpos - pos).get_norm());
			if (!sphere_test_comp(li_data.lpos, li_data.start, test_dir, radius*radius, t_val)) continue; // thick lines
			if (li_data.curr && sobjs != NULL && p2p_dist_sq(pos, li_data.lpos) >= (p2p_dist_sq(li_data.start, li_data.lpos) +
				max(0.0f, (li_data.curr->get_radius() - obj.obj->get_radius())))) continue;
		}
		fobj         = obj.obj;
		li_data.dist = dist_sq;
		if (sobjs != NULL) sobjs->push_back(obj.obj);
		if (li_data.first_only) break;
	}
}


free_obj *line_intersect_free_objects(line_int_data &li_data, int obj_types, unsigned align, bool align_only) {

	float const dmag(li_data.dir.mag());
	assert(dmag > 0.0);
	li_data.end    = li_data.start + li_data.dir*(li_data.length/dmag);
	free_obj *fobj = NULL;

	if (obj_types & OBJ_TYPE_FREE) { // limit to only ships (possibly of one alignment) if we're not looking for projectiles
		bool const find_projs((obj_types & OBJ_TYPE_PROJ) != 0), find_ships((obj_types & OBJ_TYPE_SHIP) != 0);
		if (align_only) assert(find_ships && align < NUM_ALIGNMENT);
		float const urm(find_projs ? uobj_rmax : urm_ship);
		vector<cached_obj> const &objs(find_projs ? c_uobjs : (align_only ? ships[align] : all_ships));
		line_intersect_fo_vector(li_data, objs, fobj, urm, find_ships);
	}
	if (obj_types & OBJ_TYPE_STAT) { // first only? don't know which one is first
		line_intersect_fo_vector(li_data, stat_objs, fobj, urm_static, 0); // slow if there are too many sobjs
	}
	if (fobj != NULL) li_data.dist = sqrt(li_data.dist);
	return fobj;
}


uobject *line_intersect_objects(line_int_data &li_data, free_obj *&fobj, int obj_types) {

	assert(obj_types & OBJ_TYPE_ALL);
	assert(li_data.length > 0.0);
	uobject *sobj = NULL;
	float f_dist(0.0), s_dist(0.0);
	fobj = NULL;
	
	if (obj_types & OBJ_TYPE_FOBJ) { // could allow searches for ships of some alignment through input arguments
		fobj   = line_intersect_free_objects(li_data, obj_types, 0, 0);
		f_dist = li_data.dist;
	}
	if (obj_types & OBJ_TYPE_UOBJ) { // check for universe collisions
		sobj = line_intersect_universe(li_data.start, li_data.dir, li_data.length, li_data.line_radius, s_dist);
		if (sobj != NULL && li_data.sobjs != NULL) li_data.sobjs->push_back(sobj);
	}
	if (fobj != NULL && sobj != NULL) {
		if (f_dist < s_dist) sobj = NULL; else fobj = NULL;
	}
	li_data.dist = 0.0;
	if (fobj != NULL) li_data.dist = f_dist;
	if (sobj != NULL) li_data.dist = s_dist;
	if (sobj != NULL) return sobj;
	return fobj;
}


// *********************** QUERY EXECUTORS *******************************



void apply_one_exp(query_data &qdata, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	if (cobj.obj == qdata.ptr) return; // don't apply an object's explosion to itself
	free_obj *obj(cobj.obj);
	assert(obj != NULL);
	if ((qdata.eflags & EXP_FLAGS_NO_FFIRE) && qdata.parent != NULL && qdata.parent->is_related(obj)) return;
	if ((qdata.eflags & EXP_FLAGS_NO_PART)  && (cobj.flags & OBJ_FLAGS_PART)) return;
	free_obj::intersect_params ip;
	ip.calc_dscale = 1;
	if (!obj->sphere_int_obj(qdata.pos, qdata.radius, ip)) return; // no detailed intersection
	vector3d dist(cobj.pos, qdata.pos);
	float rtot(cobj.radius);
	
	if (qdata.wclass >= 0) { // kind of a hack - add the weapon radius to the target radius if we know the weapon class
		assert((unsigned)qdata.wclass < us_weapons.size());
		rtot += us_weapons[(unsigned)qdata.wclass].radius;
	}
	float const dmag(dist.mag()), damage(qdata.damage*ip.dscale*calc_damage_scale(dmag, rtot, qdata.radius));
	if (damage <= 0.0) return;
	vector3d velocity(dist);
	if (dmag > TOLERANCE) velocity *= (0.00025/dmag);

	// If a ship A fires a rocket R, the rocket strikes another ship B and R explodes, and the blast radius damages ship C, then how would we get the pointer to ship A out of this?
	// If a ship A fires a rocket R, the rocket strikes another ship B and causes B to explode, and the blast radius damages ship C, then which ship (A or B) gets credit for damaging C?
	// OK, so lets say no one is responsible for the explosion, but the exploding object's parent (if is has one) is responsible for the damage.
	obj->collision(qdata.pos, velocity, 0.025*damage, 0.0, NULL, EXP_COLL_ELASTIC); // applies a force but does little damage
	obj->damage(damage, DAMAGE_EXP, qdata.pos, qdata.parent, qdata.wclass);
}


void apply_one_light(query_data &qdata, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	assert(cobj.obj != NULL);
	if (!cobj.obj->sphere_int_obj(qdata.pos, qdata.radius)) return; // no detailed intersection (optional test)

	if (cobj.flags & OBJ_FLAGS_SHIP) { // need a more exact intersection test (for light-emitting ships)
		float const dist(p2p_dist(cobj.pos, qdata.pos));

		if (dist > 2.0*cobj.radius) { // test for shadows
			free_obj *fobj(NULL);
			vector3d const delta((cobj.pos - qdata.pos)/dist); // start dir len cur ignore first_only check_parent [exp_r]
			line_int_data li_data(qdata.pos, delta, (dist - cobj.radius), qdata.parent, cobj.obj, 1, 0); // check_parent?
			if (line_intersect_objects(li_data, fobj, OBJ_TYPE_LGU) && fobj != qdata.parent) return;
		}
	}
	cobj.obj->add_light(qdata.index);
}


void gen_lightning_wrays(query_data &qdata, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	float const width(0.02*qdata.radius);
	assert(width > 0.0 && qdata.radius > 0.0 && cobj.radius > 0.0);
	assert(cobj.obj && qdata.parent);
	point p1(qdata.pos), p2(cobj.pos);
	vadd_rand(p1, qdata.dist*qdata.radius);
	vadd_rand(p2, 0.5*cobj.radius);
	if (!qdata.parent->is_related(cobj.obj)) add_lightning_wray(width, p1, p2);
}


void check_for_inc_proj(query_data &qdata, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	assert(cobj.obj != NULL && (cobj.flags & OBJ_FLAGS_PROJ));
	if (qdata.align != ALIGN_NEUTRAL && qdata.align != ALIGN_PIRATE &&
		cobj.obj->get_align() == (unsigned)qdata.align) return; // don't shoot a friendly shot (use is_enemy()?)
	vector3d const pdir(cobj.pos, qdata.pos);
	if (dot_product(pdir, cobj.obj->get_velocity()) >= 0.0) return; // projectile is not moving or is moving away
	float const damage(cobj.obj->damage_done());

	if (damage > qdata.damage) { // return the projectile with the highest damage (0 damage is ignored)
		qdata.ipos   = cobj.pos;
		qdata.radius = pdir.mag();
		qdata.damage = damage;
		qdata.fobj   = cobj.obj;
	}
}


bool update_min_d(closeness_data &qdata, unsigned ix) { // ships and decoys

	cached_obj const &cobj((*qdata.objs)[ix]);
	if (fabs(qdata.pos.x - cobj.pos.x) > qdata.dmin)      return 0;
	float const dist_sq(p2p_dist_sq(qdata.pos, cobj.pos));
	if (dist_sq >= qdata.dmin*qdata.dmin || dist_sq <= qdata.min_dist_sq) return 1;
	assert(cobj.flags & (OBJ_FLAGS_SHIP | OBJ_FLAGS_DECY));
	free_obj *obj(cobj.obj);
	free_obj const *const qq(qdata.questioner);
	assert(obj != NULL && qq != NULL);
	if (obj == qq)                                return 1; // no self query?
	if (obj->not_a_target() || (!qdata.friendly && !qq->target_valid(obj))) return 1; // don't attack friendlies or invalid objects
	if (qdata.req_dock && !obj->can_dock_with())  return 1; // not a dock
	if (obj->is_invisible())                      return 1; // can't see (cloaked ship)
	if (qdata.req_shields && !obj->has_shields()) return 1;

	if (cobj.flags & OBJ_FLAGS_DECY) { // decoy
		if (qq->get_root_parent() != NULL && qq->get_root_parent() == obj->get_root_parent()) return 1; // don't shoot a friendly

		if (qq->get_target() != obj && (qq->get_parent() == NULL || qq->get_parent()->get_target() != obj)) { // not a target
			unsigned const align1(qq->get_align()), align2(obj->get_align());
			if (align1 == align2 && align1 != ALIGN_NEUTRAL && align1 != ALIGN_PIRATE) return 1; // same team
		}
	}
	float dist(sqrt(dist_sq)), dscale(1.0);
	if (cobj.flags & OBJ_FLAGS_DECY) dscale *= 0.5; // higher priority
	
	if (cobj.flags & OBJ_FLAGS_SHIP) {
		if (qdata.q_dir != zero_vector && dist > cobj.radius) { // prefer ships in front of the questioner
			dscale *= (1.0 - min(0.5, 4.0*cobj.radius/dist)*dot_product(qdata.q_dir, (cobj.pos - qdata.pos))/dist);
			if (dist*dscale >= qdata.dmin*qdata.dscale) return 1;
		}
		if (obj->offense() == 0.0 || !obj->has_weapons()) { // non-offensive ships have low priority
			if (obj->offense() == 0.0) dscale *= 4.0;
			if (!obj->has_weapons())   dscale *= 4.0;
		}
		if (obj->disabled()) dscale *= 2.0; // disabled ships have lower priority
		if (dist*dscale >= qdata.dmin*qdata.dscale) return 1;
	}
	assert(dist > 0.0 && dscale > 0.0);
	vector3d const dir((cobj.pos - qdata.pos).get_norm());
	line_int_data li_data(qdata.pos, dir, dist, qq, NULL, 1, 0);
	free_obj *fobj;
	uobject const *coll_obj(line_intersect_objects(li_data, fobj, OBJ_TYPE_SOBJ));
	
	if (coll_obj != NULL && coll_obj != obj) { // can't see around a stellar object (star, planet, moon)
		//return 1;
		if (qdata.init_dmin > 0.0 && dist > 0.5*qdata.init_dmin) return 1; // less sensor range
		dscale *= 4.0;
		if (dist*dscale >= qdata.dmin*qdata.dscale) return 1; // prefer not to attack that target
	}
	qdata.dmin    = dist;//li_data.dist;
	qdata.dscale  = dscale;
	qdata.closest = obj;
	return 1;
}


bool get_all_close(all_query_data &qdata, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	if (fabs(qdata.pos.x - cobj.pos.x) > qdata.max_search_dist) return 0;

	if (dist_less_than(qdata.pos, cobj.pos, qdata.dmax_rscale*cobj.radius) && cobj.obj != qdata.questioner) {
		qdata.results.push_back(cobj.obj);
	}
	return 1;
}


void test_coll_query(query_data &qdata, unsigned ix) {

	qdata.index      = ix + 1;
	qdata.exit_query = 1;
}


// **************************** QUERY ITERATORS **************************


typedef void (*obj_query)     (query_data     &, unsigned);
typedef bool (*closest_query) (closeness_data &, unsigned);
typedef bool (*all_query)     (all_query_data &, unsigned);

unsigned const BAD_QUERY_FLAGS = (OBJ_FLAGS_BAD_ | OBJ_FLAGS_NEW_);


bool query_func_wrap(query_data &qdata, obj_query query_func, unsigned bad_flags, unsigned ix) {

	cached_obj const &cobj((*qdata.objs)[ix]);
	if (cobj.flags & (BAD_QUERY_FLAGS | bad_flags)) return 1;
	if (fabs(qdata.pos.x - cobj.pos.x) > (qdata.radius + qdata.urm)) return 0;
	float const rsum(qdata.radius + cobj.radius);
	if (fabs(qdata.pos.y - cobj.pos.y) > rsum) return 1;
	if (!dist_less_than(cobj.pos, qdata.pos, rsum)) return 1;
	query_func(qdata, ix);
	return 1;
}


inline bool query_func_wrap(closeness_data &qdata, closest_query query_func, unsigned bad_flags, unsigned ix) {

	return query_func(qdata, ix); // not all params used
}

inline bool query_func_wrap(all_query_data &qdata, all_query query_func, unsigned bad_flags, unsigned ix) {

	return query_func(qdata, ix); // not all params used
}


template<typename data_t, typename query> void find_close_objects(data_t &qdata, query query_func, unsigned bad_flags=0) {

	assert(qdata.objs != NULL);
	if (qdata.objs->empty()) return;
	unsigned const start(binary_search_pos(*(qdata.objs), qdata.pos)), nobjs(qdata.objs->size());
	assert(start <= nobjs);

	// search to left and right until dx > dmin
	for (int i = start; i >= 0 && i < (int)nobjs && !qdata.exit_query; --i) {
		if (!query_func_wrap(qdata, query_func, bad_flags, i)) break;
	}
	for (int i = start+1; i < (int)nobjs && !qdata.exit_query; ++i) {
		if (!query_func_wrap(qdata, query_func, bad_flags, i)) break;
	}
}


// ************************** QUERY DRIVERS ****************************


void apply_explosion(point const &pos, float radius, float damage, unsigned eflags, int wclass, uobject *ptr, free_obj const *parent) {

	assert(radius != 0.0 && damage >= 0.0);
	if (damage == 0.0) return;
	query_data qdata(&c_uobjs, pos, radius, uobj_rmax);
	qdata.damage = damage;
	qdata.eflags = eflags;
	qdata.wclass = wclass;
	qdata.ptr    = ptr;
	qdata.parent = parent;
	find_close_objects(qdata, apply_one_exp, (OBJ_FLAGS_NCOL | OBJ_FLAGS_NEXD | OBJ_FLAGS_NEW_));
}


void add_br_light(unsigned index, point const &pos, float radius, free_obj const *const parent) { // is parent necessary?

	assert((EXPLOSION_LIGHT + (int)NUM_EXP_LIGHTS) <= MAX_GL_LIGHT);
	if (!EXPLODE_LIGHTING || NUM_EXP_LIGHTS == 0) return;
	query_data qdata(&c_uobjs, pos, radius, uobj_rmax);
	qdata.index  = index;
	qdata.parent = parent;
	find_close_objects(qdata, apply_one_light, OBJ_FLAGS_NOLT);
}


void gen_lightning_from(point const &pos, float radius, float dist, free_obj const *src) {

	query_data qdata(&all_ships, pos, radius, urm_ship);
	qdata.parent = src;
	qdata.dist   = dist;
	find_close_objects(qdata, gen_lightning_wrays);
}


free_obj const *check_for_incoming_proj(point const &pos, int align, float dist) {

	if (urm_proj == 0.0) return NULL; // no projectiles
	query_data qdata(&coll_proj, pos, dist, urm_proj);
	qdata.align = align;
	find_close_objects(qdata, check_for_inc_proj, OBJ_FLAGS_NCOL);
	return ((qdata.radius < dist) ? qdata.fobj : NULL);
}


free_obj *free_obj::get_closest_ship(point const &pos, float min_dist, float max_dist, bool enemy, bool attack_all,
									 bool req_shields, bool decoy_tricked, bool dir_pref) const
{
	if (min_dist >= max_dist) return NULL;
	float dmin(max_dist);
	float const min_dist_sq(max(TOLERANCE, min_dist*min_dist));
	vector3d const q_dir(dir_pref ? get_dir() : zero_vector);
	closeness_data cdata(NULL, pos, dmin, min_dist_sq, this, req_shields, 0, !enemy);
	cdata.q_dir = q_dir;

	if (decoy_tricked && !decoys.empty()) {
		cdata.objs = &decoys;
		find_close_objects(cdata, update_min_d);
	}
	if (attack_all) {
		assert(enemy);
		cdata.objs = &all_ships;
		find_close_objects(cdata, update_min_d);
	}
	else {
		unsigned const alignment(get_align());
		assert(alignment < NUM_ALIGNMENT);
		bool testset[NUM_ALIGNMENT] = {0};

		switch (alignment) {
		case ALIGN_NEUTRAL:
			return NULL; // no enemies
		case ALIGN_GOV:
			if (!enemy) testset[alignment] = 1; // friend
			break;
		case ALIGN_PIRATE:
			if (!enemy) break; // no friends
			for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) testset[i] = 1; // all enemies
			break;
		case ALIGN_PLAYER:
			if (enemy && !player_enemy) break; // no enemies
		default: // ALIGN_RED, ALIGN_BLUE, etc.
			if (enemy) {
				for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {
					if (TEAM_ALIGNED(i) && i != alignment) testset[i] = 1;
				}
			}
			else {
				testset[alignment] = 1;
			}
		}
		for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {
			if (!testset[i]) continue;
			cdata.objs = &ships[i];
			find_close_objects(cdata, update_min_d);
		}
	}
	return cdata.closest;
}


free_obj *u_ship::get_closest_dock(float max_dist) const {

	closeness_data cdata(&ships[alignment], pos, max_dist, radius*radius, this, 0, 1, 1);
	find_close_objects(cdata, update_min_d);
	return cdata.closest;
}


unsigned check_for_obj_coll(point const &pos, float radius) { // test sobjs?

	query_data qdata(&c_uobjs, pos, radius, uobj_rmax);
	qdata.index = 0;
	find_close_objects(qdata, test_coll_query);
	return qdata.index;
}


void get_all_close_objects(all_query_data &qdata) {

	qdata.results.resize(0);
	find_close_objects(qdata, get_all_close);
}



