varying vec4 vertex;

void main()
{
	setup_texgen0();
	setup_texgen1();
	setup_texgen2();
	vertex = gl_Vertex;
	gl_Position = ftransform();
	set_fog_coord();
} 
