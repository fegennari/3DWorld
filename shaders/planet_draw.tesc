// define the number of CPs in the output patch
layout (vertices = 3) out;

uniform vec3 gEyeWorldPos;

// attributes of the input CPs
in vec3 WorldPos_CS[];
in vec2 TexCoord_CS[];
in vec3 Normal_CS[];

// attributes of the output CPs
out vec3 WorldPos_ES[];
out vec2 TexCoord_ES[];
out vec3 Normal_ES[];

float get_tess_level(float d1, float d2) {
	float dav = 0.5*(d1 + d2);
	if      (dav <= 2.0) {return 10.0;}
	else if (dav <= 5.0) {return 7.0;}
	else                 {return 3.0;}
}

void main() {
	// Set the control points of the output patch
	TexCoord_ES[gl_InvocationID] = TexCoord_CS[gl_InvocationID];
	Normal_ES[gl_InvocationID] = Normal_CS[gl_InvocationID];
	WorldPos_ES[gl_InvocationID] = WorldPos_CS[gl_InvocationID];

	// Calculate the distance from the camera to the three control points
	float EyeToVertexDistance0 = distance(gEyeWorldPos, WorldPos_ES[0]);
	float EyeToVertexDistance1 = distance(gEyeWorldPos, WorldPos_ES[1]);
	float EyeToVertexDistance2 = distance(gEyeWorldPos, WorldPos_ES[2]);

	// Calculate the tessellation levels
	gl_TessLevelOuter[0] = get_tess_level(EyeToVertexDistance1, EyeToVertexDistance2);
	gl_TessLevelOuter[1] = get_tess_level(EyeToVertexDistance2, EyeToVertexDistance0);
	gl_TessLevelOuter[2] = get_tess_level(EyeToVertexDistance0, EyeToVertexDistance1);
	gl_TessLevelInner[0] = gl_TessLevelOuter[2];
}                                                                                        