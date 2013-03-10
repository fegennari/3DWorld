void main()
{
	setup_texgen0();
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
} 
