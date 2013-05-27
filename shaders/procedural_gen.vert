varying vec3 vpos, normal; // world space
varying vec3 eye_norm;

void main()
{
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	normal   = normalize(gl_Normal);
	eye_norm = normalize(gl_NormalMatrix * normal);
	vpos     = gl_Vertex.xyz;
	setup_indir_lighting(vpos, normal);
} 
