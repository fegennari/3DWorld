uniform vec4 color = vec4(1.0);

void main()
{
	gl_Position   = fg_ftransform();
	gl_FrontColor = color;
}
