varying vec4 epos;
varying vec3 eye_norm;
varying vec2 tex_coord;

#ifdef USE_BUMP_MAP
#ifdef USE_TANGENT_VECTOR
attribute vec4 tangent;
varying vec4 tangent_v;

void setup_tbn() {
	// Building the matrix Eye Space -> Tangent Space
	tangent_v = vec4(normalize(gl_NormalMatrix * tangent.xyz), tangent.w);
}

#else
void setup_tbn() {} // do nothing
#endif
#endif
