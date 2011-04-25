uniform float height = 1.0;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 vertex = gl_Vertex;
	
	float motion_val = clamp((1.5 - 2.0*gl_TexCoord[0].s), 0.0, 1.0); // 1.0 for top vertex, 0.0 for bottom vertices
	// Note: grass motion amplitude should depend on dot(wind, gl_Normal), but the normal is incorrect
	float delta = get_wind_delta(vertex.xyz, gl_Color.g) * height * motion_val;
		
	// apply x/y delta but maintain the existing height
	vec3 v = normalize(vec3(delta*wind_x, delta*wind_y, height)) * height;
	v.z   -= height;
	vertex.xyz += v;
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec3 n = normalize(normal);
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_pt_light_comp(normal, epos, 0);
	if (enable_light1) color += add_pt_light_comp(normal, epos, 1);
	if (enable_light2) color += add_pt_light_comp(n, epos, 2);
	if (enable_light3) color += add_pt_light_comp(n, epos, 3);
	if (enable_light4) color += add_pt_light_comp(n, epos, 4);
	if (enable_light5) color += add_pt_light_comp(n, epos, 5);
	if (enable_light6) color += add_pt_light_comp(n, epos, 6);
	if (enable_light7) color += add_pt_light_comp(n, epos, 7);
	gl_FrontColor   = color;
	gl_FogFragCoord = length(epos.xyz);
} 
