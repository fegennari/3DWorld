varying vec4 epos;

void main()
{
	gl_Position = fg_ftransform();
	epos = gl_ModelViewMatrix * fg_Vertex;
	set_fog_coord(fg_Vertex);
} 
