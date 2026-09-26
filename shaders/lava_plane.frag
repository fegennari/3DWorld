uniform sampler2D lava_tex, noise_tex;
uniform float lifetime = 0.0;

in vec2 tc;

void main() {
	float density    = min(1.0, 8.0*lifetime) + 0.1;
	float sample_val = texture(noise_tex, 0.13*tc).r;
	if (density < sample_val) discard;
	vec4 texel   = texture(lava_tex, (tc + vec2(1.0)*0.3*lifetime)); // shift over time
	fg_FragColor = gl_Color*texel;
}
