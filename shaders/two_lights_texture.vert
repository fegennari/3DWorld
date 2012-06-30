void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec3 epos       = (gl_ModelViewMatrix * gl_Vertex).xyz;
	gl_FogFragCoord = length(epos); // set standard fog coord
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
} 
