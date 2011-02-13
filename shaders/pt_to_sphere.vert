void main()
{
	//gl_Position = ftransform();
	//scale = gl_NormalMatrix * gl_Normal; // eye space scaling
	gl_Position = gl_Vertex;
	vec3 scale = gl_Normal;
	gl_TexCoord[0].stq = scale; // FIXME: hack
	gl_FrontColor = gl_Color;
} 
