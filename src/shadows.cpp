// 3D World - OpenGL CS184 Computer Graphics Project - Mesh Shadow Generation code
// by Frank Gennari
// 3/6/06

#include "3DWorld.h"
#include "mesh.h"


extern bool combined_gu;
extern float sun_rot, moon_rot, sun_theta, moon_theta, zmin, zmax, FAR_CLIP;
extern point sun_pos, moon_pos, mesh_origin;
extern vector3d up_norm;


void update_sun_shadows() {

	// seems like this has to be first, but questionable
	if (!combined_gu && fabs(sun_rot/PI - 1.0) < 0.6 && fabs(sun_rot/PI - 1.0) > 0.4) {calc_visibility(MOON_SHADOW);}
	calc_visibility(SUN_SHADOW);
}

void fix_sun_moon_rot(float &angle) {angle = fix_angle(angle);}

void update_sun_and_moon() {

	float const radius(0.6f*(FAR_CLIP+X_SCENE_SIZE));
	fix_sun_moon_rot(sun_rot);
	fix_sun_moon_rot(moon_rot);
	light_factor = fabs(sun_rot/PI - 1.0);
	moon_pos     = mesh_origin + rtp_to_xyz(radius, moon_theta, moon_rot);
	sun_pos      = mesh_origin + rtp_to_xyz(radius, sun_theta,  sun_rot);
	point const light_pos((light_factor >= 0.5 || combined_gu) ? sun_pos : moon_pos);
	up_norm      = light_pos.get_norm();
}


bool light_valid(unsigned light_sources, int l, point &lpos) {

	if (!(light_sources & (1 << l))) return 0;
	if ((l == LIGHT_SUN && light_factor < 0.4) || (l == LIGHT_MOON && light_factor > 0.6)) return 0;
	if (!get_light_pos(lpos, l) || lpos.z < zmin) return 0;
	return 1;
}

bool light_valid_and_enabled(int l, point &lpos) {
	return (light_valid(0xFF, l, lpos) && is_light_enabled(l));
}
bool light_valid_and_enabled(int l) {point lpos; return light_valid_and_enabled(l, lpos);}

