varying vec3 vertex, epos;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	vertex = gl_Vertex.xyz;
	epos   = (gl_ModelViewMatrix * gl_Vertex).xyz;
	gl_Position = ftransform();
	set_fog_coord();
} 
