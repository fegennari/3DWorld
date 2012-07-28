void main()
{
#ifdef GEN_QUAD_TEX_COORDS
	set_tc0_from_vert_id();
#else
	gl_TexCoord[0] = gl_MultiTexCoord0;
#endif
	gl_Position = ftransform();
	calc_leaf_lighting();
} 
