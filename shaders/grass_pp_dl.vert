uniform float height = 1.0;

varying vec3 dlpos, normal, eye_norm; // world space

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0;
	vec3 gwdelta    = get_grass_wind_delta(gl_Vertex.xyz, height, gl_MultiTexCoord0.s);
	eye_norm        = length(gl_Normal) * (gl_NormalMatrix * normalize(normalize(gl_Normal) + gwdelta/height)); // eye space (not normalized)
	vec4 vertex     = gl_Vertex + vec4(gwdelta, 0.0);
	vec4 epos       = gl_ModelViewMatrix  * vertex;
	gl_Position     = gl_ProjectionMatrix * epos;
	normal          = normalize(gl_Normal);
	dlpos           = vertex.xyz;
	gl_FrontColor   = gl_Color;
	gl_FogFragCoord = length(epos.xyz);
} 
