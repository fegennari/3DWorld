uniform float animation_time = 0.0;
uniform float size_scale = 1.0;

void apply_vertex_animation(inout vec4 vertex, in vec2 tc) {
	vertex.y += (0.01/size_scale)*abs(sin(150.0*animation_time));
}
