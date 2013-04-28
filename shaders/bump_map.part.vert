varying vec4 epos;
varying vec3 eye_norm;

#ifdef USE_BUMP_MAP
attribute vec4 tangent;

varying vec2 bump_tc;
varying vec3 tangent_v, binorm_v;
varying vec3 ts_pos; // tangent space pos

void setup_tbn() {
	// Building the matrix Eye Space -> Tangent Space
	tangent_v  = normalize(gl_NormalMatrix * tangent.xyz);
	binorm_v   = cross(eye_norm, tangent_v);
	tangent_v *= tangent.w;

	ts_pos.x  = dot(epos.xyz, tangent_v);
	ts_pos.y  = dot(epos.xyz, binorm_v);
	ts_pos.z  = dot(epos.xyz, eye_norm);

	bump_tc = gl_TexCoord[0].st;
}
#endif
