attribute float point_size;

void main()
{
	gl_Position      = gl_Vertex;
	gl_FrontColor    = gl_Color;
	gl_TexCoord[7].s = point_size;
} 
