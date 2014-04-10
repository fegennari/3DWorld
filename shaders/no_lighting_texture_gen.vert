void main()
{
	setup_texgen_st();
	gl_Position   = fg_ftransform();
	gl_FrontColor = fg_Color;
} 
