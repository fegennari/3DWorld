varying vec4 epos;

void main()
{
	epos          = gl_ModelViewMatrix * gl_Vertex;
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
}
