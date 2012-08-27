uniform vec3 offset = vec3(0,0,0);

varying vec4 pos, color;

void main()
{
	pos = gl_Vertex + vec4(offset, 1.0);
	color = gl_Color;
	gl_Position = ftransform();
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz);
} 
