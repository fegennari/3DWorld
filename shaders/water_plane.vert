varying vec3 normal;
varying vec4 epos, proj_pos;

void main()
{
	setup_texgen0();
	gl_TexCoord[1] = gl_MultiTexCoord0;
	normal = gl_NormalMatrix * gl_Normal; // not normalized
	epos   = gl_ModelViewMatrix * gl_Vertex;
	proj_pos = ftransform();
	gl_Position = proj_pos;
	gl_FogFragCoord = length(epos.xyz);
}
