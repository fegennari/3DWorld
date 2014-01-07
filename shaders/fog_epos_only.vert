varying vec4 epos;

void main()
{
	gl_Position = ftransform();
	epos = gl_ModelViewMatrix * gl_Vertex;
	set_fog_coord(gl_Vertex);
} 
