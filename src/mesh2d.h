// 3D World
// by Frank Gennari
// mesh2d class definition
// 3/2/14
#ifndef _MESH2D_H_
#define _MESH2D_H_

#include "3DWorld.h"


struct mesh2d {

	float *pmap; // perturbation map
	bool *rmap;  // render map
	float *emap; // expand map
	point *ptsh; // point shift map
	unsigned size;

	unsigned get_index(unsigned s, unsigned t) const {assert(s < size && t <= size); return (s*(size+1) + t);}

public:
	float expand;

	mesh2d() : pmap(NULL), rmap(NULL), emap(NULL), ptsh(NULL), size(0), expand(0.0) {}
	//~mesh2d() {clear();}
	void clear();
	unsigned get_num()     const {assert(size > 0); return size*(size+1);} // square, with an extra row
	unsigned choose_rand() const {return (rand() % get_num());}
	void set_size(unsigned sz);
	template<typename T> void alloc_ptr(T *&p, T const val);
	void alloc_pmap();
	void alloc_rmap();
	void alloc_emap();
	void alloc_ptsh();
	void reset_pmap();
	void add_random(float mag, float min_mag, float max_mag, unsigned skipval=0);
	void mult_by(float val);
	void unset_rand_rmap(unsigned num_remove);
	void set_rand_expand(float mag, unsigned num_exp);
	void set_rand_translate(point const &tp, unsigned num_trans);
	void set_val(unsigned s, unsigned t, float val)      {assert(pmap); pmap[get_index(s, t)] = val;}
	float get_val(unsigned s, unsigned t) const          {assert(pmap); return pmap[get_index(s, t)];}
	void  set_rm(unsigned s, unsigned t, bool val)       {assert(rmap); rmap[get_index(s, t)] = val;}
	bool  get_rm(unsigned s, unsigned t) const           {assert(rmap); return rmap[get_index(s, t)];}
	void  set_em(unsigned s, unsigned t, float val)      {assert(emap); emap[get_index(s, t)] = val;}
	float get_em(unsigned s, unsigned t) const           {assert(emap); return emap[get_index(s, t)];}
	void  set_pt(unsigned s, unsigned t, point const &p) {assert(ptsh); ptsh[get_index(s, t)] = p;}
	point get_pt(unsigned s, unsigned t) const           {assert(ptsh); return ptsh[get_index(s, t)];}
	unsigned get_size() const {return size;}
	void draw_perturbed_sphere(point const &pos, float radius, int ndiv, bool tex_coord) const;
};


#endif // _MESH2D_H_

