#ifdef ENABLE_SHADOWS
uniform vec3 sun_pos;
uniform float sun_radius;
uniform vec4 shadow_casters[8];
#endif

float calc_shadow_atten(in vec3 world_space_pos)
{
	float atten = 1.0;
#ifdef ENABLE_SHADOWS
	for (int i = 0; i < num_shadow_casters; ++i) {
		atten *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, shadow_casters[i].xyz, shadow_casters[i].w);
	}
#endif
	return atten;
}

