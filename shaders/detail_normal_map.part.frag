uniform vec2 detail_normal_tex_scale = vec2(8.0);
uniform sampler2D detail_normal_tex;

in vec2 tc;

vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos, in vec3 normal, in float bump_scale) {
	vec3 tan  = normalize(cross(fg_ModelViewMatrix[0].xyz, normal));
	mat3 TBN  = transpose(mat3(cross(normal, tan), -tan, normalize(normal))); // world space {X, -Y, Z} for normal in +Z
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	vec3 nmap = texture(detail_normal_tex, detail_normal_tex_scale*tc).rgb; // scaled detail texture
#ifdef BLEND_DIST_DETAIL_NMAP // mix in a lower frequency normal map sample to break up tiling artifacts
	vec3 nmap2= texture(detail_normal_tex, detail_normal_tex_scale*tc/6.0).rgb;
	nmap = mix(nmap, nmap2, clamp((length(eye_pos) - 0.8), 0.0, 1.0));
#endif
	return normalize(mix(vec3(0,0,1), (2.0*nmap - 1.0), bump_scale));
}
