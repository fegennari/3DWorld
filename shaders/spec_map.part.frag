#ifdef USE_SPEC_MAP
uniform sampler2D spec_map;
// tex_coord comes from bump_map.part.frag

vec3 get_spec_color() {
	return texture2D(spec_map, tex_coord).rgb;
}
#endif
