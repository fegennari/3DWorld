varying vec3 normal; // in eye space

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position    = ftransform();
	gl_FrontColor  = gl_Color;
	normal    = normalize(gl_NormalMatrix * gl_Normal);
	vec4 epos = gl_ModelViewMatrix * gl_Vertex;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
} 
