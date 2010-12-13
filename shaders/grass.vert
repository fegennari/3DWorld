uniform float height = 1.0;

vec4 add_pt_light_comp(in vec3 normal, in vec4 epos, in int i) {
	vec4 color = add_light_comp(normal, i);
	float dist = distance(gl_LightSource[i].position, epos);
	float atten = 1.0 / (gl_LightSource[i].constantAttenuation +
						 gl_LightSource[i].linearAttenuation * dist +
						 gl_LightSource[i].quadraticAttenuation * dist * dist);
	return color * atten;
}

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 vertex = gl_Vertex;
	
	float motion_val = clamp((1.5 - 2.0*gl_TexCoord[0].s), 0.0, 1.0); // 1.0 for top vertex, 0.0 for bottom vertices
	// Note: grass motion amplitude should depend on dot(wind, gl_Normal), but the normal is incorrect
	float delta = get_wind_delta(vertex.xyz) * height * motion_val;
		
	// apply x/y delta but maintain the existing height
	vec3 v = normalize(vec3(delta*wind_x, delta*wind_y, height)) * height;
	v.z   -= height;
	vertex.xyz += v;
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	
	// calculate lighting: L0-L1 is directional, L2-L7 is point
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	if (enable_light2) color += add_pt_light_comp(normal, epos, 2);
	if (enable_light3) color += add_pt_light_comp(normal, epos, 3);
	if (enable_light4) color += add_pt_light_comp(normal, epos, 4);
	if (enable_light5) color += add_pt_light_comp(normal, epos, 5);
	if (enable_light6) color += add_pt_light_comp(normal, epos, 6);
	if (enable_light7) color += add_pt_light_comp(normal, epos, 7);
	gl_FrontColor   = color;
	gl_FogFragCoord = length(epos.xyz);
} 
