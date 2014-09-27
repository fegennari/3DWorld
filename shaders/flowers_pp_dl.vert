out vec3 dlpos, normal; // world space
out vec3 eye_norm;
out vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	normal        = normalize(fg_Normal);
	vec3 gwdelta  = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	eye_norm      = length(fg_Normal) * (fg_NormalMatrix * normalize(normal + gwdelta/height)); // eye space (not normalized), height comes from wind.part
	vec4 vertex   = fg_Vertex + vec4(gwdelta, 0.0);
	gl_Position   = fg_ModelViewProjectionMatrix * vertex;
	dlpos         = vertex.xyz;
	gl_FrontColor = fg_Color;
} 
