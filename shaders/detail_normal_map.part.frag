uniform vec2 detail_normal_tex_scale = vec2(8.0);
uniform sampler2D detail_normal_tex;

in vec2 tc;

vec3 nmap_texture_lookup(in vec2 nm_tc) {
#ifdef USE_TILE_BLEND_NMAP
	return ByExampleProceduralNoise(nm_tc);
#else
	return texture(detail_normal_tex, nm_tc).rgb; // scaled detail texture
#endif
}

vec3 get_bump_map_normal_dnm(in float bump_scale, in vec3 eye_pos) {
	vec3 nmap  = nmap_texture_lookup(detail_normal_tex_scale*tc);
#ifdef BLEND_DIST_DETAIL_NMAP // mix in a lower frequency normal map sample to break up tiling artifacts
	vec3 nmap2 = nmap_texture_lookup(detail_normal_tex_scale*tc/6.0);
	nmap = mix(nmap, nmap2, clamp((length(eye_pos) - 0.8), 0.0, 1.0));
#endif
	return normalize(mix(vec3(0,0,1), (2.0*nmap - 1.0), bump_scale));
}
vec3 get_bump_map_normal() {
	return get_bump_map_normal_dnm(1.0, vec3(0,0,1)); // FIXME: eye_pos is not valid
}
mat3 get_tbn(in float bscale, in vec3 n) {
	vec3 tan = normalize(cross(fg_ModelViewMatrix[0].xyz, n));
	return transpose(mat3(bscale*cross(n, tan), -tan, normalize(n))); // world space {X, -Y, Z} for normal in +Z
}
void bump_map_setup(inout vec3 light_dir, inout vec3 eye_pos, in vec3 n) {
	mat3 TBN  = get_tbn(1.0, n);
	light_dir = normalize(TBN * light_dir);
	eye_pos   = TBN * eye_pos;
}
vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos, in vec3 n, in float bump_scale) {
	bump_map_setup(light_dir, eye_pos, n);
	return get_bump_map_normal_dnm(bump_scale, eye_pos);
}
vec3 apply_bump_map_for_tbn(inout vec3 light_dir, inout vec3 eye_pos, in mat3 TBN) {
	light_dir = normalize(TBN * light_dir);
	eye_pos   = TBN * eye_pos;
	return get_bump_map_normal(); // in tangent space
}
