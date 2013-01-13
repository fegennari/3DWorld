varying vec3 vertex;
varying vec4 color;

void main()
{
	vertex = gl_Vertex.xyz;
	gl_Position = ftransform();
	set_fog_coord();
	color = gl_Color;
} 
