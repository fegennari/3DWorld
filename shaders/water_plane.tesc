// define the number of CPs in the output patch
layout (vertices = 4) out;

uniform vec3 camera_pos; // in world space

// attributes of the input CPs
in vec3 vertex[];
in vec2 tc[], tc2[];

// attributes of the output CPs
out vec3 vertex_ES[];
out vec2 tc_ES[], tc2_ES[];

float get_tess_level(float d1, float d2) {
	float dav = 0.5*(d1 + d2);
	return max(1, min(48, 200.0/dav));
}

void main() {
	// Set the control points of the output patch
	vertex_ES[gl_InvocationID] = vertex[gl_InvocationID];
	tc_ES    [gl_InvocationID] = tc    [gl_InvocationID];
	tc2_ES   [gl_InvocationID] = tc2   [gl_InvocationID];

	// Calculate the distance from the camera to the three control points
	float dist0 = distance(vertex[0], camera_pos); // SW
	float dist1 = distance(vertex[3], camera_pos); // SE
	float dist2 = distance(vertex[2], camera_pos); // NE
	float dist3 = distance(vertex[1], camera_pos); // NW

	// Calculate the tessellation levels
	gl_TessLevelOuter[0] = get_tess_level(dist0, dist1); // south
	gl_TessLevelOuter[1] = get_tess_level(dist0, dist3); // west
	gl_TessLevelOuter[2] = get_tess_level(dist2, dist3); // north
	gl_TessLevelOuter[3] = get_tess_level(dist1, dist2); // east
	gl_TessLevelInner[0] = (gl_TessLevelOuter[0] + gl_TessLevelOuter[1] + gl_TessLevelOuter[2] + gl_TessLevelOuter[3])/4;
	gl_TessLevelInner[1] = gl_TessLevelInner[0];
}                                                                                        