layout (points) in;
layout (triangle_strip, max_vertices=4) out;

uniform float size = 1.0;
in vec4 color[1]; // from VS

void main()
{
	fg_Color_vf = color[0];
	vec4 pos    = fg_ModelViewMatrix * gl_in[0].gl_Position;
	vec4 pts[4];
	pts[0] = fg_ProjectionMatrix * (pos + vec4(-size,  size, 0.0, 0.0));
	pts[1] = fg_ProjectionMatrix * (pos + vec4(-size, -size, 0.0, 0.0));
	pts[2] = fg_ProjectionMatrix * (pos + vec4( size,  size, 0.0, 0.0));
	pts[3] = fg_ProjectionMatrix * (pos + vec4( size, -size, 0.0, 0.0));
	output_textured_quad(pts);
}
