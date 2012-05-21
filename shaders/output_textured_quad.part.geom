void output_textured_quad(in vec4 pts[4])
{
	gl_Position = pts[0]; gl_TexCoord[0].st = vec2(0.0, 1.0); EmitVertex();
	gl_Position = pts[1]; gl_TexCoord[0].st = vec2(0.0, 0.0); EmitVertex();
	gl_Position = pts[2]; gl_TexCoord[0].st = vec2(1.0, 1.0); EmitVertex();
	EndPrimitive();
	gl_Position = pts[1]; gl_TexCoord[0].st = vec2(0.0, 0.0); EmitVertex();
	gl_Position = pts[2]; gl_TexCoord[0].st = vec2(1.0, 1.0); EmitVertex();
	gl_Position = pts[3]; gl_TexCoord[0].st = vec2(1.0, 0.0); EmitVertex();
	EndPrimitive();
}
