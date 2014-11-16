#ifdef USE_SPEC_MAP
uniform sampler2D spec_map;
// tc comes from bump_map.part.frag

vec3 get_spec_color() {
	return texture(spec_map, tc).rgb;
}
#endif
