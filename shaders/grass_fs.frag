uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec2 tc;
varying vec3 normal, norm_normal;
varying vec4 epos;

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;

	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos_smap_light0(norm_normal, epos);
	if (enable_light1) color += add_light_comp_pos_smap_light1(norm_normal, epos);
	if (enable_light2) color += add_pt_light_comp(normal, epos, 2);
	if (enable_light3) color += add_pt_light_comp(normal, epos, 3);
	if (enable_light4) color += add_pt_light_comp(normal, epos, 4);
	if (enable_light5) color += add_pt_light_comp(normal, epos, 5);
	if (enable_light6) color += add_pt_light_comp(normal, epos, 6);
	if (enable_light7) color += add_pt_light_comp(normal, epos, 7);
	color *= texel;
#ifndef NO_FOG
	color = apply_fog(color);
#endif
	gl_FragColor = color;
}
