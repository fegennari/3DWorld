layout (points) in;
layout (triangle_strip, max_vertices=4) out;

uniform vec2 xlate = vec2(0.0);

// from VS
in vec4 color_vs[1], vertex_vs[1];
in vec3 normal_vs[1];
in float size_vs[1];
in int dmin_vs[1];

out vec4 vertex, epos;
out vec3 eye_norm;
out vec2 tc;

void gen_vertex(in vec4 v, in vec3 en, in float ts, in float tt) {
	fg_Color_vf = color_vs [0]; // shared between vertices
	tc          = vec2(ts, tt);
	eye_norm    = en;
	vertex      = v;
	epos        = fg_ModelViewMatrix * (v + vec4(xlate, 0.0, 0.0));
	gl_Position = fg_ProjectionMatrix * epos;
	EmitVertex();
}

void main()
{
	vec3 normal = normal_vs[0];
	vec3 en     = fg_NormalMatrix * normal;
	vec4 pos    = vertex_vs[0];
	vec3 va     = vec3(0.0); // orthogonal vectors
	va[dmin_vs[0]] = 1.0;
	vec4 v2 = vec4(size_vs[0]*normalize(cross(normal, va.xyz)), 0.0);
	vec4 v1 = vec4(size_vs[0]*normalize(cross(normal, v2.xyz)), 0.0);
	gen_vertex((pos - v1 + v2), en, 0.0, 1.0);
	gen_vertex((pos - v1 - v2), en, 0.0, 0.0);
	gen_vertex((pos + v1 + v2), en, 1.0, 1.0);
	gen_vertex((pos + v1 - v2), en, 1.0, 0.0);
	EndPrimitive();
}
