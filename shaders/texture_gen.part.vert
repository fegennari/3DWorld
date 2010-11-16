void setup_texgen()
{
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	gl_TexCoord[0].s = dot(ecPosition, gl_EyePlaneS[0]);
	gl_TexCoord[0].t = dot(ecPosition, gl_EyePlaneT[0]);
	gl_TexCoord[0].p = 0.0;
	gl_TexCoord[0].q = 1.0;
}
