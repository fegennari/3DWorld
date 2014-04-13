layout (points) in;
layout (triangle_strip, max_vertices=3) out;

#ifdef SIZE_FROM_ATTRIB
in float size_val[1]; // from VS
#else
uniform float size = 1.0;
#endif

varying out vec2 tc;

void main()
{
#ifdef SIZE_FROM_ATTRIB
	float size = size_val[0];
#endif
	gl_FrontColor = gl_in[0].gl_FrontColor; // all colors are the same
	vec4 pos = fg_ModelViewMatrix * gl_in[0].gl_Position;
	
	gl_Position = fg_ProjectionMatrix * (pos + vec4(-2.0*size, -size, 0.0, 0.0));
	tc = vec2(-0.5, 0.0);
	EmitVertex();
	
	gl_Position = fg_ProjectionMatrix * (pos + vec4(0.0, 3.0*size, 0.0, 0.0));
	tc = vec2(0.5, 2.0);
	EmitVertex();
	
	gl_Position = fg_ProjectionMatrix * (pos + vec4(2.0*size, -size, 0.0, 0.0));
	tc = vec2(1.5, 0.0);
	EmitVertex();
}
