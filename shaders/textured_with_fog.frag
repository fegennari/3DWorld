uniform sampler2D tex0;
uniform float min_alpha = 0.0;
in vec2 tc;

void main()
{
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;
	fg_FragColor = gl_Color*texel;
#ifndef NO_FOG
	fg_FragColor = apply_fog(fg_FragColor);
#endif
#ifdef ENABLE_ALPHA_TO_COVERAGE
	// see https://medium.com/@bgolus/anti-aliased-alpha-test-the-esoteric-alpha-to-coverage-8b177335ae4f
	fg_FragColor.a = (fg_FragColor.a - min_alpha) / max(0.5*fwidth(texel.a), 0.0001) + 0.5;
#endif
}
