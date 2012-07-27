uniform int tc_start_ix = 0;

void main()
{
#ifdef GEN_QUAD_TEX_COORDS
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	gl_TexCoord[0].st = vec2(vec4(0,1,1,0)[tc_table_ix], vec4(0,0,1,1)[tc_table_ix]);
#else
	gl_TexCoord[0] = gl_MultiTexCoord0;
#endif
	gl_Position = ftransform();
	calc_leaf_lighting();
} 
