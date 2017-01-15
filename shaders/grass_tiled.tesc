// define the number of CPs in the output patch
layout (vertices = 3) out;

uniform mat4 fg_ViewMatrix;

// attributes of the input CPs
in vec3 vertex_from_vs[];
in vec2 tc[];
in vec4 fg_Color_vf[];

// attributes of the output CPs
out vec3 vertex_ES[];
out vec2 tc_ES[];
out vec4 fg_Color_vf_ES[];

void main() {
	// Set the control points of the output patch
	vertex_ES[gl_InvocationID] = vertex_from_vs[gl_InvocationID];
	tc_ES    [gl_InvocationID] = tc[gl_InvocationID];
	fg_Color_vf_ES[gl_InvocationID] = fg_Color_vf[gl_InvocationID];

	// Calculate the distance from the camera to the control point
	vec4 epos  = fg_ModelViewMatrix * vec4(vertex_from_vs[2], 1.0); // pick a single ref vertex, since they're very close together
	float dist = length(epos.xyz);

	// Calculate the tessellation levels
	float level = max(1, min(8, 3.0/dist));
	gl_TessLevelOuter[0] = level;
	gl_TessLevelOuter[1] = level;
	gl_TessLevelOuter[2] = 1;
	gl_TessLevelInner[0] = 2*level/3+1;
}
