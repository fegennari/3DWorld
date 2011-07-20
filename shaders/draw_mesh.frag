uniform sampler2D tex0, tex1;
varying vec3 eye, dlpos, normal; // world space
varying vec3 eye_norm;

void main()
{
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);

	vec4 lit_color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) lit_color += add_light_comp_pos_smap(eye_norm, vec4(eye,1), vec4(dlpos,1), 0);
	if (enable_light1) lit_color += add_light_comp_pos_smap(eye_norm, vec4(eye,1), vec4(dlpos,1), 1);
	lit_color = clamp(lit_color, 0.0, 1.0);

	if (enable_dlights) lit_color.rgb += add_dlights(dlpos, normalize(normal), eye); // dynamic lighting
	vec4 color = vec4((texel0.rgb * texel1.rgb * lit_color.rgb), (texel0.a * texel1.a * gl_Color.a));
	gl_FragColor = apply_fog(color); // add fog
}
