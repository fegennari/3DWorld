varying vec2 tc;

void main()
{
	tc          = gl_MultiTexCoord0;
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec3 epos       = (gl_ModelViewMatrix * gl_Vertex).xyz;
	gl_FogFragCoord = length(epos); // set standard fog coord
	vec4 color = vec4(0,0,0,1);
	if (enable_light0) color += add_light_comp0(normal);
	if (enable_light1) color += add_light_comp1(normal);
	gl_FrontColor = color;
}
