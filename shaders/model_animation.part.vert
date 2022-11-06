in ivec4 bone_ids;
in vec4 bone_weights;

uniform float animation_time  = 0.0;
uniform float animation_scale = 1.0;

// Note: color is used for debugging
void apply_vertex_animation(inout vec4 vertex, inout vec3 normal, inout vec4 color, in vec2 tc) {
	// TODO: use bone_ids and bone_weights
}
