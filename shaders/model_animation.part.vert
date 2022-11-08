layout(location = 4) in uvec4 bone_ids;
layout(location = 5) in vec4 bone_weights;

uniform float animation_time = 0.0;

const int MAX_MODEL_BONES = 200;
uniform mat4 bones[MAX_MODEL_BONES];

// Note: color is used for debugging
void apply_vertex_animation(inout vec4 vertex, inout vec3 normal, inout vec4 color, in vec2 tc) {
	//color.rgb = bone_weights.xyz;
	mat4 bone_transform = bones[bone_ids[0]] * bone_weights[0];
    bone_transform     += bones[bone_ids[1]] * bone_weights[1];
    bone_transform     += bones[bone_ids[2]] * bone_weights[2];
    bone_transform     += bones[bone_ids[3]] * bone_weights[3];
    vertex = bone_transform * vertex;
    normal = (bone_transform * vec4(normal, 1.0)).xyz;
}
