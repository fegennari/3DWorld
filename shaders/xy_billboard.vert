varying float world_space_zval;

void main()
{
	set_tc0_from_vert_id();
	world_space_zval = gl_Vertex.z;
	gl_Position     = ftransform();
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz); // set standard fog coord
	vec3 color      = vec3(0.0);
	if (enable_light0) color += add_light_comp(vec3(0,0,1), 0).rgb;
	if (enable_light1) color += add_light_comp(vec3(0,0,1), 1).rgb;
	gl_FrontColor = vec4(color, gl_Color.a);
} 
