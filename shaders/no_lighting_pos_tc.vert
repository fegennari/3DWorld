out vec2 tc;
out vec3 vertex;

void main()
{
	tc            = fg_TexCoord;
	vertex        = fg_Vertex.xyz;
	gl_Position   = fg_ftransform();
	gl_FrontColor = fg_Color;
} 
