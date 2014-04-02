varying vec3 vertex;
varying vec4 color;

void main()
{
	vertex = fg_Vertex.xyz;
	gl_Position = fg_ftransform();
	set_fog_coord(fg_Vertex);
	color = fg_Color;
} 
