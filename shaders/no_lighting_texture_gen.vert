void main()
{
	setup_texgen0();
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
} 
