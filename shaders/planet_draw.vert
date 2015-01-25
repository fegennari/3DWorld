uniform mat4 fg_ViewMatrix;
uniform float obj_radius;

out vec3 normal, world_space_pos, vertex;
out vec2 tc;

void main()
{
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
	tc     = fg_TexCoord;
	normal    = normalize(fg_NormalMatrix * fg_Normal);
	vec4 epos = fg_ModelViewMatrix * vec4(vertex, 1.0);
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position = fg_ProjectionMatrix * epos;
}

