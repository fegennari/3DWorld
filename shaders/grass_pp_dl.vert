uniform float snow_cov_amt = 0.0;

out vec3 dlpos, normal; // world space
out vec3 eye_norm;
out vec2 tc;

void main()
{
	tc           = get_grass_tc();
	normal       = normalize(fg_Normal);
	vec3 gwdelta = get_grass_wind_delta(fg_Vertex.xyz, tc.s);
	eye_norm     = normalize(fg_NormalMatrix * normalize(normal + gwdelta/height)); // eye space, height comes from wind.part
	vec4 vertex  = fg_Vertex + vec4(gwdelta, 0.0);
	gl_Position  = fg_ModelViewProjectionMatrix * vertex;
	dlpos        = vertex.xyz;
	fg_Color_vf  = mix(fg_Color, vec4(0.9, 0.9, 1.0, 1.0), snow_cov_amt);
} 
