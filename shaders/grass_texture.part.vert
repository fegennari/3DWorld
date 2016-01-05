vec2 get_grass_tc() {
	//return fg_TexCoord;
	int tix = gl_VertexID % 3;
	if (tix == 0) {return vec2(0.9, 0.1);} else if (tix == 1) {return vec2(0.9, 0.9);} else {return vec2(0.1, 0.5);} // 0.1 border around grass blade texture
} 
