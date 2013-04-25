uniform sampler2D tex0;

varying vec3 dlpos, normal; // world space
varying vec3 eye_norm;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 epos  = gl_ModelViewMatrix * vec4(dlpos, 1.0);

	// calculate lighting: L0-L1 is directional
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0 ) color += add_light_comp_pos_smap_light0(eye_norm, epos);
	if (enable_light1 ) color += add_light_comp_pos_smap_light1(eye_norm, epos);
	if (enable_dlights) color.rgb += gl_Color.rgb * add_dlights(dlpos, normal, gl_ModelViewMatrixInverse[3].xyz, vec3(1,1,1)); // dynamic lighting
	vec4 fin_color = vec4(texel.rgb * color.rgb, gl_Color.a);
#ifndef NO_FOG
	fin_color = apply_fog(fin_color);
#endif
	gl_FragColor = fin_color;
}
