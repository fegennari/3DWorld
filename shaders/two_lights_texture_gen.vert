void main()
{
	setup_texgen0();
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz); // set standard fog coord
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp0(normal);
	if (enable_light1) color += add_light_comp1(normal);
	gl_FrontColor = color;
} 
