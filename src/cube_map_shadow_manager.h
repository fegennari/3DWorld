// 3D World
// by Frank Gennari
// Cube Map Shadow Manager Classes
// 3/28/17
#pragma once

#include "3DWorld.h"

struct cube_map_lix_t {
	int ixs[6]; // one per cube face, -1 is disabled
	cube_map_lix_t() {for (unsigned i = 0; i < 6; ++i) {ixs[i] = -1;}}
	void add_cube_face_lights(point const &pos, float radius, colorRGBA const &color, float near_clip);
};

class cube_map_shadow_manager {

	typedef map<unsigned, cube_map_lix_t> obj_to_light_map_t;
	obj_to_light_map_t obj_to_light_map;
	vector<unsigned> light_free_list;

public:
	void remove_light(int light_id);
	void remove_lights(cube_map_lix_t const &lix);
	void remove_obj_light(unsigned obj_id);
	unsigned alloc_light();
	cube_map_lix_t add_obj(unsigned obj_id, bool add_light);
	void sync_light_pos(unsigned obj_id, point const &obj_pos) const;
};

