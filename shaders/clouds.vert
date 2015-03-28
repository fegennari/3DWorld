uniform vec3 camera_pos;
out vec3 vertex;
out vec4 color;

void main()
{
	vertex = fg_Vertex.xyz;
	gl_Position = fg_ftransform();
	vec3 clipped_vert; // unused
	gl_FogFragCoord = get_water_fog_coord(fg_Vertex, camera_pos, clipped_vert);
	color = fg_Color;
} 
