uniform vec2 detail_normal_tex_scale = vec2(8.0);
uniform sampler2D detail_normal_tex;

varying vec2 tc;

vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos, in vec3 normal, in float bump_scale) {
	vec3 tan  = normalize(cross(fg_ModelViewMatrix[0].xyz, normal));
	mat3 TBN  = transpose(mat3(tan, cross(normal, tan), normalize(normal))); // world space {Y, X, Z} for normal in +Z
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return normalize(mix(vec3(0,0,1), (2.0*texture2D(detail_normal_tex, detail_normal_tex_scale*tc).rgb - 1.0), bump_scale)); // scaled detail texture
}
