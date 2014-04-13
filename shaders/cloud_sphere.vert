void main()
{
	setup_texgen_st();
	gl_Position = fg_ftransform();
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec4 color  = vec4(0.0, 0.0, 0.0, fg_Color.a); // get alpha directly from fg_Color
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	if (enable_light4) color.rgb += add_pt_light_comp(normal, (fg_ModelViewMatrix * fg_Vertex), 4).rgb;
	gl_FrontColor = min(vec4(1.0), color);
} 
