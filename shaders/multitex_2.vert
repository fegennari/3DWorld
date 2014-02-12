varying vec2 tc;

void main()
{
	tc            = gl_MultiTexCoord0;
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
}
