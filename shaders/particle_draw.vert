in float point_size;

out float size_val; // to GS
out vec4 color; // to GS

void main()
{
	gl_Position = fg_Vertex;
	color       = fg_Color;
	size_val    = point_size;
} 
