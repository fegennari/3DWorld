//in vec3 vpos, normal; // world space
//in vec3 eye_norm; // comes from bump_map.frag

#ifdef USE_BUMP_MAP
uniform float bump_tex_scale = 1.0;

vec3 get_bump_map_normal() {
	return normalize(lookup_triplanar_texture(bump_tex_scale*vpos, normalize(normal), bump_map, bump_map, bump_map).xyz * 2.0 - 1.0);
	//return lookup_triplanar_texture_bump(bump_tex_scale*vpos, normalize(normal), bump_map);
}
vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
	vec3 binormal = normalize(cross(vec3(1,0,0), eye_norm));
	vec3 tangent  = normalize(cross(eye_norm, binormal));
	mat3 TBN  = transpose(mat3(tangent, binormal, eye_norm));
	light_dir = normalize(TBN * light_dir);
	eye_pos   = TBN * eye_pos;
	return get_bump_map_normal();
}
#endif // USE_BUMP_MAP

