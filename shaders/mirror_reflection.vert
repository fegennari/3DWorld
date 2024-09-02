out vec4 proj_pos;

void main() {
	proj_pos = fg_ftransform();
	gl_Position = proj_pos;
}
