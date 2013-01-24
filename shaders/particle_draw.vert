void main()
{
	gl_Position      = gl_Vertex;
	gl_FrontColor    = gl_Color;
	gl_TexCoord[7].s = gl_Normal.x; // radius encoded into normal which is placed in TC7.s
} 
