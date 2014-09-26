uniform float height_scale = 1.0;
out vec3 normal;
out vec4 epos;

void main()
{
	float hval   = gen_cloud_alpha(fg_Vertex.xyz);
	float height = height_scale * hval;
	vec4 vertex  = fg_Vertex + vec4(height*fg_Normal, 0.0);
	normal       = normalize(fg_NormalMatrix * fg_Normal);
	epos         = fg_ModelViewMatrix * vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_FrontColor= fg_Color;
} 
