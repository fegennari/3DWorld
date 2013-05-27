varying vec4 epos;
varying vec3 dlpos;
varying vec3 normal; // world space
varying vec3 eye_norm;

void main()
{
	setup_texgen0();
#ifdef HAVE_DETAIL_TEXTURE
	setup_texgen1();
#endif
	gl_Position = ftransform();
	eye_norm = gl_NormalMatrix * gl_Normal; // eye space, not normalized
	normal = normalize(gl_Normal);
	dlpos  = gl_Vertex.xyz;
	epos   = gl_ModelViewMatrix * gl_Vertex;
	gl_FrontColor   = gl_Color;
} 
