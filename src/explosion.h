// 3D World - OpenGL CS184 Computer Graphics Project - Explosion/blast radius header
// by Frank Gennari
// 10/15/05

#ifndef _EXPLOSION_H_
#define _EXPLOSION_H_


#include "3DWorld.h"

// light[EXPLOSION_LIGHT] to light[EXPLOSION_LIGHT + NUM_EXP_LIGHTS - 1] must be available (currently 2-7)
// Note: Lights 0 and 1 are for star and galaxy ambient; lights 6 and 7 are for engine lights
int      const EXPLOSION_LIGHT    = 2;
unsigned const NUM_EXP_LIGHTS     = MAX_SHADER_LIGHTS - EXPLOSION_LIGHT; // 6
int      const ENGINE_START_LIGHT = 6; // Note: if engine lights are enabled, they will overwrite the lowest influence dynamic lights
int      const ENGINE_DEF_LIGHT   = 7;
int      const BLAST_TIME         = 8; // default
unsigned const EXP_FLAGS_NO_FFIRE = 0x01;
unsigned const EXP_FLAGS_SHIP     = 0x02;
unsigned const EXP_FLAGS_NO_PART  = 0x04;


enum {ETYPE_NONE=0, ETYPE_FIRE, ETYPE_NUCLEAR, ETYPE_ENERGY, ETYPE_ATOMIC, ETYPE_PLASMA, ETYPE_EMP, ETYPE_STARB,
	  ETYPE_FUSION, ETYPE_EBURST, ETYPE_ESTEAL, ETYPE_ANIM_FIRE, ETYPE_SIEGE, NUM_ETYPES};


class free_obj; // forward references
class uobject;


struct exp_type_params {

	float duration;
	colorRGBA c1, c2;

	exp_type_params(float dur, colorRGBA const &c1_, colorRGBA const &c2_) : duration(dur), c1(c1_), c2(c2_) {}
};


struct blastr { // size = 118 (120)

	int time, st_time, type, src;
	float size, cur_size, damage;
	point pos;
	vector3d dir, up_vector;
	colorRGBA color1, color2, cur_color;
	free_obj const *parent;
	bool one_frame_only, one_frame_seen;

	blastr() {}
	blastr(int tm, int ty, int sr, float sz, float dam, point const &p, vector3d const &d,
		colorRGBA const &c1, colorRGBA const &c2, free_obj const *const pa=NULL, bool ofo=0)
		: time(tm), st_time(tm), type(ty), src(sr), size(sz), cur_size(sz), damage(dam), pos(p), dir(d.get_norm()),
		up_vector(plus_y), color1(c1), color2(c2), cur_color(c1), parent(pa), one_frame_only(ofo), one_frame_seen(0) {}
	void check_pointers();
	void update();
	void add_as_dynamic_light() const;
	bool next_frame(unsigned i);
	void process() const;
	float priority() const {return size*(time/st_time);} // times start off large and decrement to 0
};


struct explosion {

	unsigned flags;
	int time, wclass;
	float radius, intensity;
	point pos;
	uobject *source;
	free_obj const *parent;

	explosion() : flags(0), time(0), wclass(-1), radius(0.0), intensity(0.0), pos(all_zeros), source(NULL), parent(NULL) {}
	explosion(point const &pos_, float radius_, float intensity_, unsigned flags_, int wclass_, uobject *src=NULL, free_obj const *parent_=NULL)
		: flags(flags_), time(0), wclass(wclass_), radius(radius_), intensity(intensity_), pos(pos_), source(src), parent(parent_) {}
	void check_pointers();
};


inline float calc_damage_scale(float dist, float radius, float bradius) {

  return ((dist < radius) ? 1.0 : min(1.0, max(0.1, (1.0 - (dist - radius)/(bradius + TOLERANCE)))));
}


void register_explosion(point const &pos, float radius, float damage, unsigned eflags, int wclass, uobject *src, free_obj const *parent);
void add_blastr(point const &pos, vector3d const &dir, float size, float damage, int time, int src,
				colorRGBA const &color1, colorRGBA const &color2, int type=ETYPE_FIRE, free_obj const *const parent=NULL);
void setup_point_light(point const &pos, colorRGBA const &color, float radius, unsigned gl_light, shader_t *shader);
bool setup_br_light(unsigned index, point const &pos, unsigned gl_light, shader_t *shader);
bool higher_priority(unsigned first, unsigned second);


#endif


