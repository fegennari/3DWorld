#ifndef SIZE_FROM_NORMAL
uniform float size = 1.0;
#endif

void main()
{
#ifdef SIZE_FROM_NORMAL
	float size = gl_TexCoordIn[0][7].s;
#endif

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