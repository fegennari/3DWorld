layout (triangles) in;
layout (triangle_strip, max_vertices=3) out;

in vec2 tc_in[3]; // from VS
in vec4 color[3]; // from VS

out vec2 tc;

void main()
{
	for (int i = 0; i < 3; ++i) {
		fg_Color_vf = color[i];
		gl_Position = gl_in[i].gl_Position;
		tc = tc_in[i];
		EmitVertex();
	}
}
