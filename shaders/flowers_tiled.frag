uniform float min_alpha = 0.0;
uniform sampler2D tex0;

in vec3 dlpos, normal; // world space
in vec3 eye_norm;
in vec2 tc;

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;

	// calculate lighting: L0-L1 is directional
	vec4 epos  = fg_ModelViewMatrix * vec4(dlpos, 1.0);
	vec3 color = vec3(0.0);
	if (enable_light0) color += add_light_comp_pos0(eye_norm, epos).rgb;
	if (enable_light1) color += add_light_comp_pos1(eye_norm, epos).rgb;
	//if (enable_light2) color += add_light_comp_pos(eye_norm, epos, 2).rgb * calc_light_atten(epos, 2); // lightning
	fg_FragColor = apply_fog_epos(texel*vec4(color, gl_Color.a), epos);
}
