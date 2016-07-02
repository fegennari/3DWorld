//in vec3 vpos, normal; // world space
//in vec3 eye_norm; // comes from bump_map.frag

#ifdef USE_BUMP_MAP
uniform float bump_tex_scale = 1.0;

vec3 get_bump_map_normal() {
	return normalize(lookup_triplanar_texture(bump_tex_scale*vpos, normalize(normal), bump_map, bump_map, bump_map).xyz * 2.0 - 1.0);
	//return lookup_triplanar_texture_bump(bump_tex_scale*vpos, normalize(normal), bump_map);
}
mat3 get_tbn(in float bscale, in vec3 n) {
	vec3 binormal = normalize(cross(vec3(1,0,0), n));
	vec3 tangent  = normalize(cross(n, binormal));
	return transpose(mat3(tangent, bscale*binormal, n));
}
#endif // USE_BUMP_MAP

