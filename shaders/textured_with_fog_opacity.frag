uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float opacity   = 1.0;

in vec2 tc;

void main()
{
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	fg_FragColor = gl_Color*texel;
#ifndef NO_FOG
	fg_FragColor = apply_fog(fg_FragColor);
#endif
#ifdef ENABLE_ALPHA_TO_COVERAGE
	fg_FragColor.a = (fg_FragColor.a - min_alpha) / max(0.5*fwidth(texel.a), 0.0001) + 0.5;
#endif
}
