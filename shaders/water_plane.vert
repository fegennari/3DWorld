varying vec3 normal;
varying vec4 epos, proj_pos;

void main()
{
	setup_texgen(0);
	normal = gl_NormalMatrix * gl_Normal; // not normalized
	epos   = gl_ModelViewMatrix * gl_Vertex;
	proj_pos = ftransform();
	gl_Position = proj_pos;
	set_fog();
}
