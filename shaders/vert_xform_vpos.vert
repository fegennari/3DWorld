out vec3 vpos;

void main()
{
	gl_Position = fg_ftransform();
	vpos        = fg_Vertex.xyz;
}
