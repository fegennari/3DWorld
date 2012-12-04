varying vec3 world_space_pos;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0; // needed for flare texture
	world_space_pos = gl_Vertex.xyz;
	gl_Position     = ftransform();
	gl_FrontColor   = gl_Color; // not needed?
}
