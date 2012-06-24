void main()
{
	gl_FrontColor   = gl_Color;
	gl_Position     = gl_Vertex;
	vec4 epos       = gl_ModelViewMatrix * gl_Vertex;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	gl_TexCoord[7]  = gl_MultiTexCoord7; // height is TC7.s - FIXME
} 
