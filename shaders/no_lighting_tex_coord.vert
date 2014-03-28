uniform vec3 xlate = vec3(0);
uniform vec3 scale = vec3(1);

varying vec2 tc;

void main()
{
	tc            = gl_MultiTexCoord0;
	gl_Position   = gl_ModelViewProjectionMatrix * (vec4(xlate, 0.0) + (vec4(scale, 1.0) * gl_Vertex));
	gl_FrontColor = gl_Color;
} 
