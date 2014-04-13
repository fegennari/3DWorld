varying vec2 tc;

void main()
{
	tc          = get_grass_tc();
	vec3 gwdelta= get_grass_wind_delta(fg_Vertex.xyz, tc.s);
	vec3 n      = fg_NormalMatrix * normalize(normalize(fg_Normal) + gwdelta/height); // eye space (not normalized), height comes from wind.part
	vec3 normal = n*length(fg_Normal); // convert to original mag (for shadows)
	vec4 vertex = fg_Vertex + vec4(gwdelta, 0.0);
	vec4 epos   = fg_ModelViewMatrix  * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec4 color = vec4(0,0,0,1);
	if (enable_light0) color.rgb += add_light_comp_pos_smap_light0(normal, epos).rgb;
	if (enable_light1) color.rgb += add_light_comp_pos_smap_light1(normal, epos).rgb;
	if (enable_light2) color.rgb += add_pt_light_comp(n, epos, 2).rgb;
	if (enable_light3) color.rgb += add_pt_light_comp(n, epos, 3).rgb;
	if (enable_light4) color.rgb += add_pt_light_comp(n, epos, 4).rgb;
	if (enable_light5) color.rgb += add_pt_light_comp(n, epos, 5).rgb;
	if (enable_light6) color.rgb += add_pt_light_comp(n, epos, 6).rgb;
	if (enable_light7) color.rgb += add_pt_light_comp(n, epos, 7).rgb;
	gl_FrontColor   = color;
#ifndef NO_FOG
	gl_FogFragCoord = length(epos.xyz);
#endif
} 
