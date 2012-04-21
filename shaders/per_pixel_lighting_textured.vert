varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}
	vec3 n = gl_NormalMatrix * gl_Normal;
	normal = (no_normalize ? n : normalize(n));
	epos   = gl_ModelViewMatrix * gl_Vertex;
	dlpos  = gl_Vertex.xyz;
	eye    = gl_ModelViewMatrixInverse[3].xyz; // world space
	gl_Position = ftransform();
	set_fog();
}
