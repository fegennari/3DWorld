out vec4 epos, proj_pos;

void main() {
	epos        = fg_ModelViewMatrix * fg_Vertex;
	proj_pos    = fg_ProjectionMatrix * epos;
	gl_Position = proj_pos;
	fg_Color_vf = fg_Color;
}
