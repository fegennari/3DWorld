void main()
{
#ifdef GEN_QUAD_TEX_COORDS
	set_tc0_from_vert_id();
#else
	tc = fg_TexCoord;
#endif
	gl_Position = fg_ftransform();
	calc_leaf_lighting();
} 
