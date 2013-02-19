varying vec3 vertex;
varying vec4 epos;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	vertex = gl_Vertex.xyz;
	epos   = (gl_ModelViewMatrix * gl_Vertex);
	gl_Position = ftransform();
	set_fog_coord();
} 
