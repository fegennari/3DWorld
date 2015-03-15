out float world_space_zval;

void main()
{
	set_tc0_from_vert_id();
	world_space_zval = fg_Vertex.z;
	vec4 epos        = fg_ModelViewMatrix * fg_Vertex;
	gl_Position      = fg_ProjectionMatrix * epos;
	gl_FogFragCoord  = length(epos.xyz); // set standard fog coord
	vec3 color       = vec3(0.0);
	if (enable_light0) color += add_light_comp_pos0(vec3(0,0,1), epos).rgb;
	if (enable_light1) color += add_light_comp_pos1(vec3(0,0,1), epos).rgb;
	fg_Color_vf = vec4(color, fg_Color.a);
} 
