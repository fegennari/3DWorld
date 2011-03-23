varying vec3 normal;

void main()
{
	normal = normalize(gl_NormalMatrix * gl_Normal);
	gl_Position = ftransform();
	set_fog();
}
