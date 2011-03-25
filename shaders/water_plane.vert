varying vec3 normal;
varying vec4 epos;

void main()
{
	setup_texgen(0);
	normal = gl_NormalMatrix * gl_Normal; // not normalized
	epos   = gl_ModelViewMatrix * gl_Vertex;
	gl_Position = ftransform();
	set_fog();
}
