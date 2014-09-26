uniform vec4 texgen_s, texgen_t;
out vec2 tc;

void setup_texgen_st() {
	tc = vec2(dot(fg_Vertex, texgen_s), dot(fg_Vertex, texgen_t));
}

uniform int tc_start_ix = 0;

void set_tc0_from_vert_id() {
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	tc = vec2(vec4(0,1,1,0)[tc_table_ix], vec4(0,0,1,1)[tc_table_ix]);
} 
