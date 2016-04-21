layout (points) in;
layout (triangle_strip, max_vertices=4) out;

uniform vec3 camera_pos;
uniform vec3 up_vector;

in vec4 vertex_vs[1];
in vec4 color_vs[1];
in vec2 delta_vs[1];

out vec4 eye_space_pos;
#ifdef TREE_BRANCHES
out vec4 world_space_pos;
#endif
out vec2 tc;

void do_vertex(in vec4 pos, in vec3 delta, in float ts, in float tt) {
	fg_Color_vf  = color_vs[0];
	vec4 out_pos = pos + vec4(delta, 0.0);
#ifdef TREE_BRANCHES
	world_space_pos = out_pos;
#endif
	eye_space_pos = fg_ModelViewMatrix * out_pos;
	gl_Position   = fg_ProjectionMatrix * eye_space_pos;
	tc = vec2(ts, tt);
	EmitVertex();
}

void main() {

	vec4 pos    = vertex_vs[0];
#ifndef TREE_BRANCHES // leaves
	pos.xyz    += 0.5*delta_vs[0].x*normalize(camera_pos - pos.xyz);
#endif
	vec3 vdir   = camera_pos - pos.xyz; // z
	vec3 v1     = normalize(cross(vdir, up_vector))*delta_vs[0].x; // what if colinear?
#ifdef TREE_BRANCHES
	vec3 v2     = delta_vs[0].y*up_vector;
#else // leaves
	vec3 v2     = delta_vs[0].y*normalize(cross(v1, vdir));
#endif
	do_vertex(pos,  v1-v2, 1.0, 0.0);
	do_vertex(pos, -v1-v2, 0.0, 0.0);
	do_vertex(pos,  v1+v2, 1.0, 1.0);
	do_vertex(pos, -v1+v2, 0.0, 1.0);
	EndPrimitive();
}
