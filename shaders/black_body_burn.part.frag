uniform sampler2D burn_mask;
uniform float burn_offset    = 0.0;
uniform float burn_tex_scale = 1.0;

vec4 get_black_body_color(in vec4 base_color, in float burn_level) {
	float v = clamp(burn_level, 0.0, 1.0);
	return mix(base_color, vec4(clamp(2.0*(v-0.25)*(v-0.25)*(v-0.25), 0, 1), clamp(500.0*(v-0.9)*(v-0.9)*(v-0.9), 0, 1), 0.0, 1.0), v);
}

vec4 apply_black_body_burn_mask(in vec4 base_color, in vec2 tc) {
	float burn_level = texture(burn_mask, burn_tex_scale*tc).r + burn_offset;
	if (burn_level >= 1.0) discard; // optional?
	return get_black_body_color(base_color, burn_level);
}

