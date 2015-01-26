uniform float obj_radius;

out vec3 vertex;

void main() {
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
}

