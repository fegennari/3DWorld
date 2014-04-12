varying vec4 epos;
varying vec3 dlpos;
varying vec3 normal; // world space
varying vec3 eye_norm;

void main()
{
	setup_texgen_st();
	eye_norm = gl_NormalMatrix * fg_Normal; // eye space, not normalized
	normal = normalize(fg_Normal);
	dlpos  = fg_Vertex.xyz;
	epos   = gl_ModelViewMatrix * fg_Vertex;
	gl_Position   = gl_ProjectionMatrix * epos;
	gl_FrontColor = fg_Color;
} 
