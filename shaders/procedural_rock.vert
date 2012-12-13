uniform float height_scale = 1.0;
varying vec3 normal;
varying vec4 epos;

void main()
{
	float hval   = gen_cloud_alpha(gl_Vertex.xyz);
	float height = height_scale * hval;
	vec4 vertex  = gl_Vertex + vec4(height*gl_Normal, 0.0);
	normal       = normalize(gl_NormalMatrix * gl_Normal);
	epos         = gl_ModelViewMatrix * vertex;
	gl_Position     = gl_ProjectionMatrix * epos;
	gl_FrontColor   = gl_Color;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
} 
