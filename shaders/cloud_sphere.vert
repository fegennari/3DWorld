void main()
{
	setup_texgen0();
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec4 color  = vec4(0.0, 0.0, 0.0, gl_Color.a); // get alpha directly from gl_Color
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	if (enable_light4) color.rgb += add_pt_light_comp(normal, (gl_ModelViewMatrix * gl_Vertex), 4).rgb;
	gl_FrontColor = min(vec4(1.0), color);
} 
