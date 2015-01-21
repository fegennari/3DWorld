layout(triangles, equal_spacing, ccw) in;

uniform mat4 gVP;
uniform sampler2D gDisplacementMap;
uniform float gDispFactor;

in vec3 WorldPos_ES[];
in vec2 TexCoord_ES[];
in vec3 Normal_ES[];

out vec3 WorldPos_FS;
out vec2 TexCoord_FS;
out vec3 Normal_FS;

vec2 interpolate2D(vec2 v0, vec2 v1, vec2 v2) {
	return vec2(gl_TessCoord.x) * v0 + vec2(gl_TessCoord.y) * v1 + vec2(gl_TessCoord.z) * v2;
}
vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2) {
	return vec3(gl_TessCoord.x) * v0 + vec3(gl_TessCoord.y) * v1 + vec3(gl_TessCoord.z) * v2;
}

void main() {
	// Interpolate the attributes of the output vertex using the barycentric coordinates
	TexCoord_FS = interpolate2D(TexCoord_ES[0], TexCoord_ES[1], TexCoord_ES[2]);
	Normal_FS = interpolate3D(Normal_ES[0], Normal_ES[1], Normal_ES[2]);
	Normal_FS = normalize(Normal_FS);
	WorldPos_FS = interpolate3D(WorldPos_ES[0], WorldPos_ES[1], WorldPos_ES[2]);

	// Displace the vertex along the normal
	float Displacement = texture(gDisplacementMap, TexCoord_FS.xy).x;
	WorldPos_FS += Normal_FS * Displacement * gDispFactor;
	gl_Position = gVP * vec4(WorldPos_FS, 1.0);
}
