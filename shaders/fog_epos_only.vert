out vec4 epos;

void main()
{
	epos = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	set_fog_coord(fg_Vertex);
} 
