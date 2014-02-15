#ifdef ENABLE_INSTANCING
attribute vec3 xlate;
uniform float vertex_scale = 1.0;
#endif

varying float world_space_zval;

void main()
{
	set_tc0_from_vert_id();
	vec4 vertex = gl_Vertex;
#ifdef ENABLE_INSTANCING
	vertex.xyz *= vertex_scale;
	vertex.xyz += xlate;
#endif
	world_space_zval = vertex.z;
	vec4 epos        = gl_ModelViewMatrix  * vertex;
	gl_Position      = gl_ProjectionMatrix * epos;
	gl_FogFragCoord  = length(epos.xyz); // set standard fog coord
	vec3 normal      = normalize(gl_NormalMatrix * gl_Normal);
	//normal          *= normal.z/abs(normal.z); // two-sided lighting
	vec3 color       = vec3(0.0);
	if (enable_light0) {color += add_light_comp0(normal).rgb;}
	if (enable_light1) {color += add_light_comp1(normal).rgb;}
	if (enable_light2) {color += add_pt_light_comp(normal, epos, 2).rgb;}
	gl_FrontColor = vec4(color, gl_Color.a);
}
