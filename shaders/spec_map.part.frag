// Note: tc comes from bump_map.part.frag
#ifdef USE_SPEC_MAP
uniform sampler2D spec_map;
vec3 get_spec_color() {return texture(spec_map,  tc).rgb;}
#endif

#ifdef USE_GLOSS_MAP
uniform sampler2D gloss_map;
float get_gloss_val() {return texture(gloss_map, tc).r  ;} // 0.0-1.0
#endif
