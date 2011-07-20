varying vec3 eye, dlpos, normal; // world space
varying vec3 eye_norm;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	eye_norm = gl_NormalMatrix * gl_Normal; // eye space
	normal = normalize(gl_Normal);
	dlpos  = gl_Vertex.xyz;
	eye    = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	gl_FrontColor = gl_Color;
	set_fog();
} 
