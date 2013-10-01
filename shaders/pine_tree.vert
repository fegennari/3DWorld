void main()
{
	set_tc0_from_vert_id();
	gl_Position     = ftransform();
	vec4 epos       = gl_ModelViewMatrix * gl_Vertex;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	vec3 normal     = normalize(gl_NormalMatrix * gl_Normal);
	//normal         *= normal.z/abs(normal.z); // two-sided lighting
	vec3 color      = vec3(0.0);
	if (enable_light0) {color += add_light_comp(normal, 0).rgb;}
	if (enable_light1) {color += add_light_comp(normal, 1).rgb;}
	if (enable_light2) {color += add_pt_light_comp(normal, epos, 2).rgb;}
	gl_FrontColor = vec4(color, gl_Color.a);
}
