void main()
{
	setup_texgen(0);
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
} 
