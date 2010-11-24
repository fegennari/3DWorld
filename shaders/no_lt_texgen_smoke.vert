varying vec3 eye, vpos;

void main()
{
	setup_texgen(0);
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
	eye  = (inverse(gl_ModelViewMatrix) * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	vpos = gl_Vertex.xyz;
} 
