uniform float wind_time  = 0.0;
uniform float wind_mag   = 0.0;
uniform float wind_scale = 1.0; // used as 0/1 toggle
uniform float wind_add   = 0.0;
uniform float wind_freq  = 1.0;
uniform float wind_zscale= 1.0;
//uniform float tex_coord_weight; // from leaf_wind.part

void add_leaf_wind(inout vec4 vpos) {
#ifdef ENABLE_WIND
	int tc_table_ix = (gl_VertexID + tc_start_ix) & 3;
	if (tex_coord_weight == 0.0 && tc_table_ix >= 2) return; // only move the tip verts (not the end verts)
	float wmag = wind_mag * wind_scale + wind_add; // * (1.0 - fg_TexCoord.t);

	if (wmag > 0.0) { // only move the tip verts (not the end verts)
		float t   = wind_time;
		vec3 pos  = wind_freq*vpos.xyz;
		vpos.xyz += wmag*vec3((sin(1.00*t + pos.x) + sin(2.75*t + 1.5*pos.x)),
		                      (sin(1.13*t + pos.y) + sin(2.66*t + 1.5*pos.y)),
							  (sin(1.24*t + pos.z) + sin(2.52*t + 1.5*pos.z))*wind_zscale);
	}
#endif // ENABLE_WIND
}
