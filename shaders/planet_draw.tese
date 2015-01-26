layout(triangles, equal_spacing, ccw) in;

uniform float disp_scale = 1.0;
uniform sampler2D displace_map;

in vec3 vertex_ES[];
in vec3 normal_ES[];

out vec3 vertex;
out vec3 normal;

vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2) {
	return vec3(gl_TessCoord.x) * v0 + vec3(gl_TessCoord.y) * v1 + vec3(gl_TessCoord.z) * v2;
}

void main() {
	// Interpolate the attributes of the output vertex using the barycentric coordinates
	normal = normalize(interpolate3D(normal_ES[0], normal_ES[1], normal_ES[2]));
	vertex = interpolate3D(vertex_ES[0], vertex_ES[1], vertex_ES[2]);

	// Displace the vertex along the normal
	//float displacement = texture(displace_map, TexCoord_FS.xy).x;
	//vertex += normal * displacement * disp_scale;
	gl_Position = fg_ModelViewProjectionMatrix * vec4(vertex, 1.0);
}
