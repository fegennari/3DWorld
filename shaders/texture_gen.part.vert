varying vec2 tc, tc2, tc3;

// Note: probably okay to only copy the s and t components
void setup_texgen0() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	tc = vec4(dot(ecPosition, gl_EyePlaneS[0]), dot(ecPosition, gl_EyePlaneT[0]), 0.0, 1.0);
}

void setup_texgen1() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	tc2 = vec4(dot(ecPosition, gl_EyePlaneS[1]), dot(ecPosition, gl_EyePlaneT[1]), 0.0, 1.0);
}

void setup_texgen2() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	tc3 = vec4(dot(ecPosition, gl_EyePlaneS[2]), dot(ecPosition, gl_EyePlaneT[2]), 0.0, 1.0);
}

uniform int tc_start_ix = 0;

void set_tc0_from_vert_id()
{
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	tc = vec2(vec4(0,1,1,0)[tc_table_ix], vec4(0,0,1,1)[tc_table_ix]);
} 
