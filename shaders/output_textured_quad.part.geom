out vec2 tc;

void output_textured_quad(in vec4 pts[4], in vec4 color_out) { // all colors are the same
	fg_Color_vf = color_out; gl_Position = pts[0]; tc = vec2(0.0, 1.0); EmitVertex();
	fg_Color_vf = color_out; gl_Position = pts[1]; tc = vec2(0.0, 0.0); EmitVertex();
	fg_Color_vf = color_out; gl_Position = pts[2]; tc = vec2(1.0, 1.0); EmitVertex();
	fg_Color_vf = color_out; gl_Position = pts[3]; tc = vec2(1.0, 0.0); EmitVertex();
	EndPrimitive();
}
