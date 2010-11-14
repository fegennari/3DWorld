#version 120
#extension GL_EXT_geometry_shader4 : enable

uniform float size;
 
// a passthrough geometry shader for color and position
void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos = gl_ModelViewMatrix * gl_PositionIn[0];
	
	gl_Position = gl_ProjectionMatrix * (pos + vec4(-2.0*size, -size, 0.0, 0.0));
	gl_TexCoord[0].st = vec2(-0.5, 0.0);
	EmitVertex();
	
	gl_Position = gl_ProjectionMatrix * (pos + vec4(0.0, 3.0*size, 0.0, 0.0));
	gl_TexCoord[0].st = vec2(0.5, 2.0);
	EmitVertex();
	
	gl_Position = gl_ProjectionMatrix * (pos + vec4(2.0*size, -size, 0.0, 0.0));
	gl_TexCoord[0].st = vec2(1.5, 0.0);
	EmitVertex();
}