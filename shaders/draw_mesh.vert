varying vec4 epos, dlpos;
varying vec3 normal; // world space
varying vec3 eye_norm;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	eye_norm = gl_NormalMatrix * gl_Normal; // eye space
	normal = normalize(gl_Normal);
	dlpos  = gl_Vertex;
	epos   = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor = gl_Color;
	set_fog();
} 
