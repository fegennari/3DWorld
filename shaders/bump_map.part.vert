#ifdef USE_BUMP_MAP
#ifdef USE_TANGENT_VECTOR
in vec4 tangent;
out vec4 tangent_v;

void setup_tbn() {
	// Building the matrix Eye Space -> Tangent Space
	tangent_v = vec4(normalize(fg_NormalMatrix * tangent.xyz), tangent.w);
}

#else
void setup_tbn() {} // do nothing
#endif
#endif
