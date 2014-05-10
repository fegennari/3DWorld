void main()
{
	setup_texgen_st();
	gl_Position = fg_ftransform();
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec3 color  = vec3(0.0); // get alpha directly from fg_Color
	if (enable_light0) color += add_light_comp0(normal).rgb;
	if (enable_light1) color += add_light_comp1(normal).rgb;
	if (enable_light4) color += add_pt_light_comp(normal, (fg_ModelViewMatrix * fg_Vertex), 4).rgb;
	gl_FrontColor = vec4(min(vec3(1.0), color), fg_Color.a);
} 
