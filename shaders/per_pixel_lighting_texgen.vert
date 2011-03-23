varying vec3 normal;

void main()
{
	setup_texgen(0);
	normal = normalize(gl_NormalMatrix * gl_Normal);
	gl_Position = ftransform();
	set_fog();
}
