varying vec4 epos;

void main()
{
	epos = gl_ModelViewMatrix * fg_Vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	set_fog_coord(fg_Vertex);
} 
