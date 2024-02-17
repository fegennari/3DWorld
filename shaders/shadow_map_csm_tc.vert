out vec2 tc_vs;

void main() {
	gl_Position = fg_Vertex; // pass-through with no transforms
	tc_vs       = fg_TexCoord;
}
