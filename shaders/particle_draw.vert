attribute float point_size;

void main()
{
	gl_Position      = fg_Vertex;
	gl_FrontColor    = fg_Color;
	gl_TexCoord[7].s = point_size;
} 
