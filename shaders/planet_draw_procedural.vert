uniform float obj_radius;

out vec3 normal, vertex;

void main()
{
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
	normal = displace_vertex_and_get_normal(fg_Vertex.xyz, vertex);
	gl_Position = fg_ModelViewProjectionMatrix * vec4(vertex, 1.0);
}

