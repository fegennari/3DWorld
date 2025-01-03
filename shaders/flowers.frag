uniform float min_alpha    = 0.0;
uniform float snow_cov_amt = 0.0;
uniform sampler2D tex0;

in vec3 dlpos, normal; // world space
in vec3 eye_norm;
in vec2 tc;

void main()
{
	//if (tc.s <= 0.0 || tc.s >= 1.0 || tc.t <= 0.0 || tc.t >= 1.0) discard;
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;

	// calculate lighting: L0-L1 is directional
	vec4 epos  = fg_ModelViewMatrix * vec4(dlpos, 1.0);
	vec3 color = vec3(0.0);
	if (enable_light0 ) color += add_light_comp_pos_smap_light0(normalize(eye_norm), epos).rgb;
	if (enable_light1 ) color += add_light_comp_pos_smap_light1(normalize(eye_norm), epos).rgb;
	if (enable_dlights) add_dlights(color, dlpos, epos, normalize(normal), gl_Color.rgb); // dynamic lighting
	fg_FragColor = mix(texel, vec4(1.0), snow_cov_amt)*vec4(color, gl_Color.a);
#ifndef NO_FOG
	fg_FragColor = apply_fog_epos(fg_FragColor, epos);
#endif
}
