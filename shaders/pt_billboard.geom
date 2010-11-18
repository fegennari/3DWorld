#version 120
#extension GL_EXT_geometry_shader4 : enable

uniform float size = 1.0;
 
// a passthrough geometry shader for color and position
void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	
	vec4 pos = gl_ModelViewMatrix * gl_PositionIn[0];
	vec4 pts[4];
	pts[0] = gl_ProjectionMatrix * (pos + vec4(-size,  size, 0.0, 0.0));
	pts[1] = gl_ProjectionMatrix * (pos + vec4(-size, -size, 0.0, 0.0));
	pts[2] = gl_ProjectionMatrix * (pos + vec4( size,  size, 0.0, 0.0));
	pts[3] = gl_ProjectionMatrix * (pos + vec4( size, -size, 0.0, 0.0));
	
	gl_Position = pts[0]; gl_TexCoord[0].st = vec2(0.0, 1.0); EmitVertex();
	gl_Position = pts[1]; gl_TexCoord[0].st = vec2(0.0, 0.0); EmitVertex();
	gl_Position = pts[2]; gl_TexCoord[0].st = vec2(1.0, 1.0); EmitVertex();
	EndPrimitive();
	
	gl_Position = pts[1]; gl_TexCoord[0].st = vec2(0.0, 0.0); EmitVertex();
	gl_Position = pts[2]; gl_TexCoord[0].st = vec2(1.0, 1.0); EmitVertex();
	gl_Position = pts[3]; gl_TexCoord[0].st = vec2(1.0, 0.0); EmitVertex();
	EndPrimitive();
}