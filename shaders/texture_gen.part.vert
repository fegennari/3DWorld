uniform vec4 texgen_s, texgen_t;
out vec2 tc;

void setup_texgen_st() {
	tc = vec2(dot(fg_Vertex, texgen_s), dot(fg_Vertex, texgen_t));
}

uniform int tc_start_ix = 0;

void set_tc0_from_vert_id() {
	int tix = (gl_VertexID + tc_start_ix) & 3;
	if (tix == 0) {tc = vec2(0,0);} else if (tix == 1) {tc = vec2(1,0);} else if (tix == 2) {tc = vec2(1,1);} else {tc = vec2(0,1);}
} 
