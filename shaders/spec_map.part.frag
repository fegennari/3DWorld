#ifdef USE_SPEC_MAP
uniform sampler2D spec_map;

vec3 get_spec_color() {
	return texture2D(spec_map, gl_TexCoord[0].st).rgb;
}
#endif
