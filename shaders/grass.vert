varying vec2 tc;

void main()
{
	tc          = get_grass_tc();
	vec3 gwdelta= get_grass_wind_delta(gl_Vertex.xyz, tc.s);
	vec3 n      = gl_NormalMatrix * normalize(normalize(gl_Normal) + gwdelta/height); // eye space (not normalized), height comes from wind.part
	vec3 normal = n*length(gl_Normal); // convert to original mag (for shadows)
	vec4 vertex = gl_Vertex + vec4(gwdelta, 0.0);
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos_smap_light0(normal, epos);
	if (enable_light1) color += add_light_comp_pos_smap_light1(normal, epos);
	if (enable_light2) color += add_pt_light_comp(n, epos, 2);
	if (enable_light3) color += add_pt_light_comp(n, epos, 3);
	if (enable_light4) color += add_pt_light_comp(n, epos, 4);
	if (enable_light5) color += add_pt_light_comp(n, epos, 5);
	if (enable_light6) color += add_pt_light_comp(n, epos, 6);
	if (enable_light7) color += add_pt_light_comp(n, epos, 7);
	gl_FrontColor   = color;
#ifndef NO_FOG
	gl_FogFragCoord = length(epos.xyz);
#endif
} 
