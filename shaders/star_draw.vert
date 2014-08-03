uniform float radius = 1.0;

varying vec3 world_space_pos;
varying vec2 tc;

void main()
{
	tc              = fg_TexCoord; // needed for flare texture
	world_space_pos = radius*fg_Vertex.xyz;
	gl_Position     = fg_ftransform();
}
