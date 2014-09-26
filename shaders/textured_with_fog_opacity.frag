uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float opacity   = 1.0;

in vec2 tc;

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	fg_FragColor = apply_fog(gl_Color*texel);
}
