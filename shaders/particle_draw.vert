attribute float point_size;
out float size_val; // to GS

void main()
{
	gl_Position   = fg_Vertex;
	gl_FrontColor = fg_Color;
	size_val      = point_size;
} 
