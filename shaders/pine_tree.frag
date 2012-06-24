uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec3 normal; // in eye space

void main()
{
	vec4 texel  = texture2D(tex0, gl_TexCoord[0].st);
	float alpha = texel.a * gl_Color.a;
#ifndef NO_ALPHA_TEST
	if (alpha <= min_alpha) discard;
#endif
	// directional light sources with no attenuation
	vec3 lit_color = vec3(0,0,0);
	vec3 n = ((!two_sided_lighting || gl_FrontFacing) ? normalize(normal) : -normalize(normal)); // two-sided lighting
	if (enable_light0) lit_color += add_light_comp(n, 0).rgb;
	if (enable_light1) lit_color += add_light_comp(n, 1).rgb;
	vec4 color = vec4((texel.rgb * lit_color), alpha);
#ifndef NO_FOG
	color = apply_fog(color);
#endif
	gl_FragColor = color;
}
