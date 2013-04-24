uniform float height = 1.0;
varying vec2 tc;
varying vec3 normal, norm_normal;
varying vec4 epos;

void main()
{
	tc          = gl_MultiTexCoord0;
	vec3 gwdelta= get_grass_wind_delta(gl_Vertex.xyz, height, gl_MultiTexCoord0.s);
	normal      = gl_NormalMatrix * normalize(normalize(gl_Normal) + gwdelta/height); // eye space (not normalized)
	norm_normal = normal*length(gl_Normal); // convert to original mag (for shadows)
	vec4 vertex = gl_Vertex + vec4(gwdelta, 0.0);
	epos        = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	gl_FrontColor   = gl_Color;
	gl_FogFragCoord = length(epos.xyz);
} 
