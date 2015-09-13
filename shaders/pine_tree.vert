#ifdef ENABLE_INSTANCING
in vec3 xlate;
uniform float vertex_scale = 1.0;
#endif

out float world_space_zval;

void main()
{
	set_tc0_from_vert_id();
	vec4 vertex = fg_Vertex;
#ifdef ENABLE_INSTANCING
	vertex.xyz *= vertex_scale;
	vertex.xyz += xlate;
#endif
	add_leaf_wind(vertex);
	world_space_zval = vertex.z;
	vec4 epos        = fg_ModelViewMatrix  * vertex;
	gl_Position      = fg_ProjectionMatrix * epos;
	gl_FogFragCoord  = length(epos.xyz); // set standard fog coord
	vec3 normal      = normalize(fg_NormalMatrix * fg_Normal);
	//normal          *= normal.z/abs(normal.z); // two-sided lighting
	vec3 color       = vec3(0.0);
	if (enable_light0) {color += add_light_comp0(normal).rgb;}
	if (enable_light1) {color += add_light_comp1(normal).rgb;}
	if (enable_light2) {color += add_pt_light_comp(normal, epos, 2).rgb;}
	fg_Color_vf = vec4(color, fg_Color.a);
}
