uniform float height = 1.0;
varying vec2 tc;

void main()
{
	tc          = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 vertex = gl_Vertex;
	vertex.xyz += get_grass_wind_delta(vertex.xyz, height, gl_MultiTexCoord0.s);
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec3 n = normalize(normal);
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
	gl_FogFragCoord = length(epos.xyz);
} 
