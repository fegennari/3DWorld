out vec3 world_space_pos;
out vec2 tc;

void main()
{
	tc              = fg_TexCoord; // needed for flare texture
	world_space_pos = fg_Vertex.xyz;
	gl_Position     = fg_ftransform();
}
