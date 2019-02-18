//uniform float min_alpha = 0.0;
in vec2 tc;

void main() {
	vec4 texel = vec4(ByExampleProceduralNoise(tc), 1);
	//if (texel.a <= min_alpha) discard;
	fg_FragColor = gl_Color*texel;
#ifndef NO_FOG
	fg_FragColor = apply_fog(fg_FragColor);
#endif
}
