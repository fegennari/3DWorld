void main()
{
	gl_Position   = ftransform();
	gl_FrontColor = vec4(gl_Normal, 1.0); // world space (not normalized)
} 
