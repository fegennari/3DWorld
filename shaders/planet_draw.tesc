// define the number of CPs in the output patch
layout (vertices = 3) out;

uniform vec3 camera; // in world space

// attributes of the input CPs
in vec3 vertex[];
in vec3 normal[];

// attributes of the output CPs
out vec3 vertex_ES[];
out vec3 normal_ES[];

float get_tess_level(float d1, float d2) {
	float dav = 0.5*(d1 + d2);
	if      (dav <= 2.0) {return 10.0;}
	else if (dav <= 5.0) {return 7.0;}
	else                 {return 3.0;}
}

void main() {
	// Set the control points of the output patch
	normal_ES[gl_InvocationID] = normal[gl_InvocationID];
	vertex_ES[gl_InvocationID] = vertex[gl_InvocationID];

	// Calculate the distance from the camera to the three control points
	float dist0 = distance(camera, vertex_ES[0]);
	float dist1 = distance(camera, vertex_ES[1]);
	float dist2 = distance(camera, vertex_ES[2]);

	// Calculate the tessellation levels
	gl_TessLevelOuter[0] = get_tess_level(dist1, dist2);
	gl_TessLevelOuter[1] = get_tess_level(dist2, dist0);
	gl_TessLevelOuter[2] = get_tess_level(dist0, dist1);
	gl_TessLevelInner[0] = gl_TessLevelOuter[2];
}                                                                                        