varying vec3 normal;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	normal = normalize(gl_NormalMatrix * gl_Normal);
	gl_Position = ftransform();
	set_fog();
}
