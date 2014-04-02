uniform float height_scale = 1.0;
varying vec3 normal;
varying vec4 epos;

void main()
{
	float hval   = gen_cloud_alpha(fg_Vertex.xyz);
	float height = height_scale * hval;
	vec4 vertex  = fg_Vertex + vec4(height*fg_Normal, 0.0);
	normal       = normalize(gl_NormalMatrix * fg_Normal);
	epos         = gl_ModelViewMatrix * vertex;
	gl_Position  = gl_ProjectionMatrix * epos;
	gl_FrontColor= fg_Color;
} 
