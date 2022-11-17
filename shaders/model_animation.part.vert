layout(location = 4) in uvec4 bone_ids;
layout(location = 5) in vec4 bone_weights;

const int MAX_MODEL_BONES = 200; // must agree with the value used in model3d::setup_bone_transforms()
uniform mat4 bones[MAX_MODEL_BONES];

void apply_vertex_animation(inout vec4 vertex, inout vec3 normal, in vec2 tc) {
	mat4 bone_transform = bones[bone_ids[0]] * bone_weights[0];
    bone_transform     += bones[bone_ids[1]] * bone_weights[1];
    bone_transform     += bones[bone_ids[2]] * bone_weights[2];
    bone_transform     += bones[bone_ids[3]] * bone_weights[3];
    vertex = bone_transform * vertex;
    normal = inverse(transpose(mat3(bone_transform))) * normal;
}
