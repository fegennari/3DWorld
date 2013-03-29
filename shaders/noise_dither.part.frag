uniform sampler2D noise_tex;
uniform float noise_tex_size = 128.0;

void check_noise_and_maybe_discard(in float min_noise, in float max_noise) {
	float noise_val = texture2D(noise_tex, 0.5*gl_FragCoord.xy/noise_tex_size).r;
	if (noise_val < min_noise || noise_val > max_noise) discard;
}

void check_noise_and_maybe_discard_tc(in float min_noise, in float max_noise, in vec2 tc) { // Note: unused
	float noise_val = texture2D(noise_tex, gl_TexCoord[0].st).r;
	if (noise_val < min_noise || noise_val > max_noise) discard;
}
