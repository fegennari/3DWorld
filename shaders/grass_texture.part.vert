vec2 get_grass_tc() {
	//return fg_TexCoord;
	int tc_table_ix = gl_VertexID % 3;
	return vec2(vec3(0.9,0.9,0.1)[tc_table_ix], vec3(0.1,0.9,0.5)[tc_table_ix]); // 0.1 border around grass blade texture
} 
