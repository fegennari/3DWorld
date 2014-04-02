varying vec3 normal, vertex; // local object space

void main()
{
	normal       = fg_Normal;
	vertex       = fg_Vertex.xyz;
	gl_Position  = fg_ftransform();
}

