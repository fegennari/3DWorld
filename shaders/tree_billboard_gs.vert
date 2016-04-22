out vec4 vertex_vs;
out vec4 color_vs;
out vec2 delta_vs;

void main() {
	vertex_vs = fg_Vertex;
	color_vs  = fg_Color;
	delta_vs  = fg_TexCoord.st; // Note: could use vec2 delta attribute
}
