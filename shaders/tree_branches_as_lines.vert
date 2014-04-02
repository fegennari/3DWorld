uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
attribute float radius;

void main()
{
	gl_Position   = fg_Vertex;
	gl_FrontColor = fg_Color;
	// use fg_Normal?
	gl_TexCoord[7].r = radius;
} 
