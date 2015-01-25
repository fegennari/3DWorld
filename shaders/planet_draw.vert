uniform float obj_radius;

out vec3 normal, vertex;
out vec2 tc;

void main()
{
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
	tc     = fg_TexCoord;
	normal = normalize(fg_NormalMatrix * fg_Normal);
	gl_Position = fg_ModelViewProjectionMatrix * vec4(vertex, 1.0);
}

