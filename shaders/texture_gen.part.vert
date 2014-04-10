varying vec2 tc;

// Note: probably okay to only copy the s and t components
void setup_texgen0() {
	vec4 epos = gl_ModelViewMatrix * fg_Vertex;
	tc = vec2(dot(epos, gl_EyePlaneS[0]), dot(epos, gl_EyePlaneT[0]));
}

uniform int tc_start_ix = 0;

void set_tc0_from_vert_id()
{
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	tc = vec2(vec4(0,1,1,0)[tc_table_ix], vec4(0,0,1,1)[tc_table_ix]);
} 
