uniform sampler2D noise_tex;
uniform float noise_tex_size = 128.0;

void check_noise_and_maybe_discard(in float min_noise, in float max_noise) {
	vec2 tc = 2.0*(gl_FragCoord.xy + vec2(0.5, 0.5))/noise_tex_size;
	float noise_val = texture(noise_tex, tc).r;
	if (noise_val < min_noise || noise_val > max_noise) discard;
}
