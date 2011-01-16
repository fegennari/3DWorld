varying vec3 eye, dlpos, normal;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	vec3 norm = gl_NormalMatrix * gl_Normal; // eye space, not normalized
	normal = normalize(gl_Normal);
	dlpos  = gl_Vertex.xyz;
	eye    = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(norm, 0);
	if (enable_light1) color += add_light_comp(norm, 1);
	gl_FrontColor = clamp(color, 0.0, 1.0);
	set_fog();
} 
