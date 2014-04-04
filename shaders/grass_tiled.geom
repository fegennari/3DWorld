layout (triangles) in;
layout (triangle_strip, max_vertices=3) out;

in vec2 tc_in[3]; // from VS
varying out vec2 tc;

void main()
{
	gl_FrontColor = gl_in[0].gl_FrontColor; // all colors are the same

	for (int i = 0; i < 3; ++i) {
		gl_Position = gl_in[i].gl_Position;
		tc = tc_in[i];
		EmitVertex();
	}
}
