layout (points) in;
layout (triangle_strip, max_vertices=4) out;

uniform float radius_scale = 1.0;
uniform vec2 xlate;

in vec4 vertex_vs[1];
in vec4 color_vs[1];
in vec2 xz_sz_vs[1];

out vec2 tc;
out float world_space_zval;

void do_vertex(in vec4 pos, in vec3 delta, in float ts, in float tt) {
	fg_Color_vf = color_vs[0];
	vec4 vertex = pos + vec4(delta, 0.0);
	world_space_zval = vertex.z;
	vec4 epos   = fg_ModelViewMatrix * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	tc = vec2(ts, tt);
	EmitVertex();
}

void main() {
	vec4 pos = vertex_vs[0];
	vec3 dir = normalize(vec3(-(pos.y + xlate.y), (pos.x + xlate.x), 0.0)); // cross(z, pos-camera_pos)
	vec3 dx  = radius_scale * xz_sz_vs[0].x * dir;
	vec3 dz  = vec3(0.0, 0.0, xz_sz_vs[0].y); // for two top verts
	do_vertex(pos,  dx,    1.0, 1.0);
	do_vertex(pos, -dx,    0.0, 1.0);
	do_vertex(pos,  dx+dz, 1.0, 0.0);
	do_vertex(pos, -dx+dz, 0.0, 0.0);
	EndPrimitive();
}

