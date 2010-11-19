void setup_texgen(int i)
{
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	gl_TexCoord[i].s = dot(ecPosition, gl_EyePlaneS[i]);
	gl_TexCoord[i].t = dot(ecPosition, gl_EyePlaneT[i]);
	gl_TexCoord[i].p = 0.0;
	gl_TexCoord[i].q = 1.0;
}
