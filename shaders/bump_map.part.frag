varying vec4 epos;
varying vec3 eye_norm;
varying vec2 bump_tc;

#ifdef USE_BUMP_MAP
uniform sampler2D bump_map;

varying vec3 tangent_v, binorm_v;
varying vec3 ts_pos; // tangent space pos

vec3 apply_tbn(in vec3 v) {
	vec3 new_v;
	new_v.x = dot(v, tangent_v);
	new_v.y = dot(v, binorm_v);
	new_v.z = dot(v, eye_norm);
	return new_v;
}

// Note: we use tex coord 0 for the bump map (assuming it's the same tex coord for the regular texture)
vec3 apply_bump_map(inout vec3 light_dir) {
	light_dir = apply_tbn(light_dir);
	return normalize(texture2D(bump_map, bump_tc).xyz * 2.0 - 1.0);
}
#endif
