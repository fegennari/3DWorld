uniform sampler2D burn_mask;
uniform float burn_offset    = 0.0;
uniform float burn_tex_scale = 1.0;

vec4 get_black_body_color(in vec4 base_color, in float burn_level) {
	float v = clamp(burn_level, 0.0, 1.0);
	//return mix(base_color, vec4(min(1.0, 2.0*v), v, 0.5*v*v, 1.0), v);
	return mix(base_color, vec4(clamp(8.0*abs(v-0.5)*(v-0.5), 0, 1), clamp(75.0*abs(v-0.9)*(v-0.9), 0, 1), 0.0, 1.0), v);
}

vec4 apply_burn_mask(in vec4 base_color, in vec2 tc) {
	float burn_level = texture2D(burn_mask, burn_tex_scale*tc).a + burn_offset;
	if (burn_level >= 1.0) discard; // optional?
	return get_black_body_color(base_color, burn_level);
}

