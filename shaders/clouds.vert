varying vec3 vertex;
varying vec4 color;

void main()
{
	vertex = gl_Vertex.xyz;
	gl_Position = ftransform();
	set_fog_coord(gl_Vertex);
	color = gl_Color;
} 
