uniform vec4 color = vec4(1.0);

void main()
{
	gl_Position   = ftransform();
	gl_FrontColor = color;
}
