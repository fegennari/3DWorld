uniform float snow_cov_amt = 0.0;

out vec3 dlpos, normal; // world space
out vec3 eye_norm;

void main()
{
	set_tc0_from_vert_id();
	normal       = normalize(fg_Normal);
	vec3 gwdelta = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	eye_norm     = length(fg_Normal) * (fg_NormalMatrix * normalize(normal + gwdelta/height)); // eye space (not normalized), height comes from wind.part
	vec4 vertex  = fg_Vertex + vec4(gwdelta, 0.0);
	gl_Position  = fg_ModelViewProjectionMatrix * vertex;
	dlpos        = vertex.xyz;
	fg_Color_vf  = mix(fg_Color, vec4(0.9, 0.9, 1.0, 1.0), snow_cov_amt);
} 
