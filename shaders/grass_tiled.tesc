// define the number of CPs in the output patch
layout (vertices = 3) out;

uniform float min_tess_level = 1.0;
uniform float tess_lod_scale = 1.0;

// attributes of the input CPs
in vec3 vertex_from_vs[];
in vec2 tc[];
in vec4 fg_Color_vf[];
#ifdef INCLUDE_NORMALS
in vec3 normal[], eye_norm[];
#endif

// attributes of the output CPs
out vec3 vertex_ES[];
out vec2 tc_ES[];
out vec4 fg_Color_vf_ES[];
#ifdef INCLUDE_NORMALS
out vec3 normal_ES[], eye_norm_ES[];
#endif

void main() {
	// Set the control points of the output patch
	vertex_ES     [gl_InvocationID] = vertex_from_vs[gl_InvocationID];
	tc_ES         [gl_InvocationID] = tc            [gl_InvocationID];
	fg_Color_vf_ES[gl_InvocationID] = fg_Color_vf   [gl_InvocationID];
#ifdef INCLUDE_NORMALS
	normal_ES     [gl_InvocationID] = normal        [gl_InvocationID];
	eye_norm_ES   [gl_InvocationID] = eye_norm      [gl_InvocationID];
#endif

	// Calculate the distance from the camera to the control point
	vec4 epos  = fg_ModelViewMatrix * vec4(vertex_from_vs[2], 1.0); // pick a single ref vertex, since they're very close together
	float dist = length(epos.xyz);

	// Calculate the tessellation levels
	// Multiply tess level by alpha so that disabled/transparent grass blades use min_tess_level? Seems to cause level mismatch for verts: fg_Color_vf[gl_InvocationID].a
	float level = max(min_tess_level, min(8, tess_lod_scale*(2.5/dist)));
	gl_TessLevelOuter[0] = level;
	gl_TessLevelOuter[1] = level;
	gl_TessLevelOuter[2] = 1;
	gl_TessLevelInner[0] = 2*level/3+1;
}
