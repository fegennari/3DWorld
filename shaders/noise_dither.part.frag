#ifndef DISABLE_NOISE_DITHER
uniform sampler2D noise_tex;
uniform float noise_tex_size = 128.0;
#endif

void check_noise_and_maybe_discard(in float min_noise, in float max_noise) {
#ifndef DISABLE_NOISE_DITHER
	float noise_val = texture2D(noise_tex, gl_FragCoord.xy/noise_tex_size).r;
	if (noise_val < min_noise || noise_val > max_noise) discard;
#endif
}
