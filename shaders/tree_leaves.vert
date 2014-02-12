void main()
{
#ifdef GEN_QUAD_TEX_COORDS
	set_tc0_from_vert_id();
#else
	tc = gl_MultiTexCoord0.st;
#endif
	gl_Position = ftransform();
	calc_leaf_lighting();
} 
