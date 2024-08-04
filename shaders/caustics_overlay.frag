uniform float time       = 0.0;
uniform float smap_scale = 1.0;
uniform float tc_scale   = 1.0;
uniform vec4 color_scale = vec4(1.0);
uniform sampler2D caustic_tex;

in vec2 tc;
in vec4 epos;
in vec3 normal; // eye space

void main() {
	vec2 tc_s    = tc_scale*tc;
	float ntime  = 2.0*abs(fract(0.005*time) - 0.5);
	vec4 color   = color_scale*mix(texture(caustic_tex, tc_s), texture(caustic_tex, (tc_s + vec2(0.3, 0.6))), ntime);
	vec4 lpos    = fg_LightSource[0].position;
	color.rgb   *= max(0.0, dot(normalize(lpos.xyz), normal));
	if (smap_scale > 0.0) {color.rgb *= mix(1.0, get_shadow_map_weight_light0(epos, normal), smap_scale);}
	fg_FragColor = color;
}
