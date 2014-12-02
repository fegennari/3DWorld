void main()
{
	setup_texgen_st();
	gl_Position = fg_ftransform();
	fg_Color_vf = fg_Color;
} 
