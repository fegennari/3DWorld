uniform vec4 color_scale = vec4(1.0);
uniform float snow_cov_amt = 0.0;
uniform sampler2D tex0;

in vec3 dlpos, normal; // world space
in vec3 eye_norm;
in vec2 tc;

void main()
{
	// calculate lighting: L0-L1 is directional
	vec4 epos  = fg_ModelViewMatrix * vec4(dlpos, 1.0);
	vec3 color = vec3(0.0);
	if (enable_light0 ) color += add_light_comp_pos_smap_light0(eye_norm, epos).rgb;
	if (enable_light1 ) color += add_light_comp_pos_smap_light1(eye_norm, epos).rgb;
	if (enable_dlights) add_dlights(color, dlpos, normal, gl_Color.rgb); // dynamic lighting
	vec4 fcolor = color_scale*vec4(color, gl_Color.a);
#ifndef NO_GRASS_TEXTURE
	fcolor *= texture(tex0, tc);
#endif
	if (snow_cov_amt > 0.0) {
		// add in snow on top of grass, using ratio of lit color from base color to pick up lighting
		// FIXME: incorrect for specular component, which shouldn't be divided by the color
		fcolor.rgb = mix(fcolor.rgb, vec3(0.9, 0.9, 1.0)*color.rgb/max(vec3(0.01), gl_Color.rgb), snow_cov_amt);
	}
#ifndef NO_FOG
	fcolor = apply_fog_epos(fcolor, epos);
#endif
	fg_FragColor = fcolor;
}
