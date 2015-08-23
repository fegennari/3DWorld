uniform float wind_time = 0.0;
uniform float wind_mag  = 0.0;
uniform float wind_freq = 1.0;

void main() {
	set_tc0_from_vert_id();
	vec4 vpos = fg_Vertex;

#ifdef ENABLE_WIND
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;

	if (wind_mag > 0.0 && tc_table_ix < 2) { // only move the tip verts (not the end verts)
		float t   = 0.1*wind_time;
		vec3 pos  = wind_freq*vpos.xyz;
		vpos.xyz += wind_mag*vec3(sin(1.00*t + pos.x) + sin(2.75*t + 1.5*pos.x),
		                          sin(1.13*t + pos.y) + sin(2.66*t + 1.5*pos.y),
								  sin(1.24*t + pos.z) + sin(2.52*t + 1.5*pos.z));
	}
#endif
	gl_Position = fg_ModelViewProjectionMatrix * vpos;
	calc_leaf_lighting();
} 
