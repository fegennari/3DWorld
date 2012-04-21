varying vec4 epos;
varying vec3 normal; // world space

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	normal = normalize(gl_NormalMatrix * gl_Normal);
	epos   = gl_ModelViewMatrix * gl_Vertex;
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
	//gl_BackColor = gl_Color;
}
