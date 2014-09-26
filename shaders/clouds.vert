uniform vec3 camera_pos;
out vec3 vertex;
out vec4 color;

void main()
{
	vertex = fg_Vertex.xyz;
	gl_Position = fg_ftransform();
	gl_FogFragCoord = get_water_fog_coord(fg_Vertex, camera_pos);
	color = fg_Color;
} 
