#define NUM_CASCADES 4

layout(triangles, invocations = 4) in;
layout(triangle_strip, max_vertices = 3) out;

uniform mat4 light_space_matrices[NUM_CASCADES];

in vec2 tc_vs[3]; // from VS
out vec2 tc;

void main() {
    for (int i = 0; i < 3; ++i) {
        gl_Position = light_space_matrices[gl_InvocationID] * gl_in[i].gl_Position;
        gl_Layer = gl_InvocationID;
        tc = tc_vs[i];
        EmitVertex();
    }
    EndPrimitive();
}
