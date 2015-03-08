out vec3 vertex;
out vec2 tc2;

void main() {
	setup_texgen_st();
	tc2    = fg_TexCoord;
	vertex = fg_Vertex.xyz;
}
