uniform int tc_start_ix = 0;

void set_tc0_from_vert_id()
{
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	gl_TexCoord[0].st = vec2(vec4(0,1,1,0)[tc_table_ix], vec4(0,0,1,1)[tc_table_ix]);
} 
