varying vec4 epos;
varying vec3 eye_norm;
varying vec2 tex_coord;

#ifdef USE_BUMP_MAP
attribute vec4 tangent;

varying vec3 tangent_v;
varying vec3 ts_pos; // tangent space pos

void setup_tbn() {
	// Building the matrix Eye Space -> Tangent Space
	tangent_v   = normalize(gl_NormalMatrix * tangent.xyz);
	vec3 binorm = cross(eye_norm, tangent_v);
	tangent_v  *= tangent.w;
	ts_pos      = vec3(dot(epos.xyz, tangent_v), dot(epos.xyz, binorm), dot(epos.xyz, eye_norm));
}
#endif
