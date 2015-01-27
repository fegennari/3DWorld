layout(triangles, equal_spacing, ccw) in;

uniform float obj_radius;

in vec3 vertex_ES[];

out vec3 vertex;
out vec3 normal;

vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2) {
	return vec3(gl_TessCoord.x) * v0 + vec3(gl_TessCoord.y) * v1 + vec3(gl_TessCoord.z) * v2;
}

void main() {
	// Interpolate the attributes of the output vertex using the barycentric coordinates
	vec3 norm_vertex = normalize(interpolate3D(vertex_ES[0], vertex_ES[1], vertex_ES[2]));
	vertex = obj_radius*norm_vertex; // renormalize to planet radius
	normal = displace_vertex_and_get_normal(norm_vertex, vertex);
	gl_Position = fg_ModelViewProjectionMatrix * vec4(vertex, 1.0);
}
