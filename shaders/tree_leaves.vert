void main()
{
#ifdef GEN_QUAD_TEX_COORDS
	const vec4 tc_s_table = vec4(0,0,1,1);
	const vec4 tc_t_table = vec4(1,0,0,1);
	int tc_table_ix = gl_VertexID & 3;
	gl_TexCoord[0].st = vec2(tc_s_table[tc_table_ix], tc_t_table[tc_table_ix]);
#else
	gl_TexCoord[0] = gl_MultiTexCoord0;
#endif
	gl_Position = ftransform();
	calc_leaf_lighting();
} 
