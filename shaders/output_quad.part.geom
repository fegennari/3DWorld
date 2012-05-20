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


void output_quad(in vec4 pts[4], in vec3 n[4], in vec2 tc[2])
{
	gl_Position = pts[0]; normal = n[0]; gl_TexCoord[0].st = tc[0]; EmitVertex();
	gl_Position = pts[1]; normal = n[1]; gl_TexCoord[0].st = tc[1]; EmitVertex();
	gl_Position = pts[2]; normal = n[2]; gl_TexCoord[0].st = tc[2]; EmitVertex();
	EndPrimitive();
	gl_Position = pts[1]; normal = n[1]; gl_TexCoord[0].st = tc[1]; EmitVertex();
	gl_Position = pts[2]; normal = n[2]; gl_TexCoord[0].st = tc[2]; EmitVertex();
	gl_Position = pts[3]; normal = n[3]; gl_TexCoord[0].st = tc[3]; EmitVertex();
	EndPrimitive();
}