uniform float detail_tex_scale = 1.0;
uniform sampler2D tex0, tex1;

varying vec4 epos;
varying vec3 dlpos;
varying vec3 normal; // world space
varying vec3 eye_norm;

void main()
{
	vec4 lit_color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) lit_color += add_light_comp_pos_smap_light0(eye_norm, epos);
	if (enable_light1) lit_color += add_light_comp_pos_smap_light1(eye_norm, epos);
	lit_color = clamp(lit_color, 0.0, 1.0);
	if (enable_dlights) lit_color.rgb += add_dlights(dlpos, normalize(normal), gl_ModelViewMatrixInverse[3].xyz, vec3(1.0)); // dynamic lighting
	lit_color *= texture2D(tex0, gl_TexCoord[0].st);

#ifdef MULT_DETAIL_TEXTURE // for mesh
	lit_color *= texture2D(tex1, gl_TexCoord[1].st);
#endif
#ifdef ADD_DETAIL_TEXTURE // for water
	float slope_scale = clamp(50.0*(1.0 - normalize(normal).z), 0.0, 1.0); // normal.z is typically in +z
	lit_color *= vec4(1,1,1,1) + slope_scale*detail_tex_scale*texture2D(tex1, gl_TexCoord[1].st);
	lit_color  = clamp(lit_color, 0.0, 1.0);
#endif

	vec4 color = vec4(lit_color.rgb, gl_Color.a);
#ifndef NO_FOG
	color = apply_fog(color);
#endif
	gl_FragColor = color;
}
