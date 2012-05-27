uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
attribute float radius;

void main()
{
	gl_Position   = gl_Vertex;
	gl_FrontColor = gl_Color;
	// use gl_Normal?
	gl_TexCoord[7].r = radius;
} 
