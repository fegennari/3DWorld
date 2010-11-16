void main()
{
	setup_texgen();
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	gl_FrontColor = gl_Color;
} 
