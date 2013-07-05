varying vec4 epos;
varying vec3 eye_norm;
varying vec2 tex_coord;

#ifdef USE_BUMP_MAP
uniform sampler2D bump_map;

varying vec4 tangent_v;

// Note: we assume the bump map tex coords are the same as the object diffuse tex coords
vec3 apply_bump_map(inout vec3 light_dir, out vec3 eye_pos) {
	mat3 TBN  = transpose(mat3(tangent_v.xyz*tangent_v.w, cross(eye_norm, tangent_v.xyz), eye_norm));
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return normalize(texture2D(bump_map, tex_coord).xyz * 2.0 - 1.0);
}
#endif
