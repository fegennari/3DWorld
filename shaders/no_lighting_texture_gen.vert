void main()
{
	setup_texgen0();
	gl_Position = fg_ftransform();
	gl_FrontColor = fg_Color;
} 
