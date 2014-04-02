varying vec3 dlpos, normal; // world space
varying vec3 eye_norm;
varying vec2 tc;

void main()
{
	tc            = get_grass_tc();
	vec3 gwdelta  = get_grass_wind_delta(fg_Vertex.xyz, tc.s);
	eye_norm      = length(fg_Normal) * (gl_NormalMatrix * normalize(normalize(fg_Normal) + gwdelta/height)); // eye space (not normalized), height comes from wind.part
	vec4 vertex   = fg_Vertex + vec4(gwdelta, 0.0);
	gl_Position   = gl_ModelViewProjectionMatrix * vertex;
	normal        = normalize(fg_Normal);
	dlpos         = vertex.xyz;
	gl_FrontColor = fg_Color;
} 
