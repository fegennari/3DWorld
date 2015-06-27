out vec2 tc;

void main()
{
	tc          = get_grass_tc();
	vec3 gwdelta= get_grass_wind_delta(fg_Vertex.xyz, tc.s);
	vec3 normal = normalize(fg_NormalMatrix * (normalize(fg_Normal) + gwdelta/height)); // eye space, height comes from wind.part
	vec4 vertex = fg_Vertex + vec4(gwdelta, 0.0);
	vec4 epos   = fg_ModelViewMatrix  * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec3 color = vec3(0.0);
	if (enable_light0) color += add_light_comp_pos_smap_light0(normal, epos).rgb;
	if (enable_light1) color += add_light_comp_pos_smap_light1(normal, epos).rgb;
	if (enable_light2) color += add_pt_light_comp(n, epos, 2).rgb;
	fg_Color_vf = vec4(color, 1.0);
#ifndef NO_FOG
	gl_FogFragCoord = length(epos.xyz);
#endif
} 
