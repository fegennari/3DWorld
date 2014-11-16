uniform float detail_tex_scale = 1.0;
uniform sampler2D tex0, tex1;

//in vec2 tc; // comes from detail_normal_map.part.frag
//in vec4 epos; // predeclared earlier for dynamic lighting to work
in vec3 dlpos;
in vec3 normal; // world space
in vec3 eye_norm;

vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
	return apply_bump_map(light_dir, eye_pos, eye_norm, 1.0);
}

void main()
{
	vec3 lit_color = vec3(0.0);
	if (enable_light0) lit_color += add_light_comp_pos_smap_light0(eye_norm, epos).rgb;
	if (enable_light1) lit_color += add_light_comp_pos_smap_light1(eye_norm, epos).rgb;
	lit_color = clamp(lit_color, 0.0, 1.0);
	if (enable_dlights) lit_color.rgb += add_dlights(dlpos, normalize(normal), vec3(1.0)); // dynamic lighting
	lit_color *= texture(tex0, tc).rgb;

#ifdef MULT_DETAIL_TEXTURE // for mesh
	lit_color *= texture(tex1, 32.0*tc).rgb; // 32x scale
#endif
#ifdef ADD_DETAIL_TEXTURE // for water
	float slope_scale = clamp(50.0*(1.0 - normalize(normal).z), 0.0, 1.0); // normal.z is typically in +z
	lit_color *= vec3(1.0) + slope_scale*detail_tex_scale*texture(tex1, 32.0*tc).rgb; // 32x scale
	lit_color  = clamp(lit_color, 0.0, 1.0);
#endif

	fg_FragColor = vec4(lit_color, gl_Color.a);
#ifndef NO_FOG
	fg_FragColor = apply_fog_epos(fg_FragColor, epos);
#endif
}
