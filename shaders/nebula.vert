varying vec3 normal, vertex; // local object space

void main()
{
	normal       = gl_Normal;
	vertex       = gl_Vertex.xyz;
	gl_Position  = ftransform();
}

