out vec2 tc;

void main() {
	gl_Position = fg_ftransform();
	tc          = fg_TexCoord;
}
