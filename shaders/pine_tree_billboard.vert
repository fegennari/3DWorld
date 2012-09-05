void main()
{
	set_tc0_from_vert_id();
	gl_Position     = ftransform();
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz); // set standard fog coord
	vec3 nz     = (gl_NormalMatrix * vec3(0.0, 0.0, gl_Normal.z));
	vec3 normal = normalize(nz + vec3(0.0, 0.0, sqrt(1.0 - nz.x*nz.x - nz.y*nz.y - nz.z*nz.z)));
	vec3 color  = vec3(0,0,0);
	if (enable_light0) color += add_light_comp(normal, 0).rgb;
	if (enable_light1) color += add_light_comp(normal, 1).rgb;
	gl_FrontColor = vec4(color, gl_Color.a);
} 
