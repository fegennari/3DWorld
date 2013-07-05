varying vec4 epos;
varying vec3 eye_norm;
varying vec2 tex_coord;

#ifdef USE_BUMP_MAP
uniform sampler2D bump_map;

varying vec3 tangent_v;
varying vec3 ts_pos; // tangent space pos
varying float tangent_w;

vec3 apply_tbn(in vec3 v) {
	return vec3(dot(v, tangent_v), dot(v, tangent_w*cross(eye_norm, tangent_v)), dot(v, eye_norm));
}

// Note: we assume the bump map tex coords are the same as the object diffuse tex coords
vec3 apply_bump_map(inout vec3 light_dir, out vec3 eye_pos) {
	light_dir = apply_tbn(light_dir);
	eye_pos   = ts_pos; // convert to tangent space
	return normalize(texture2D(bump_map, tex_coord).xyz * 2.0 - 1.0);
}
#endif
