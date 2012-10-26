uniform vec3 cloud_offset = vec3(0,0,0);

varying vec4 pos, color;

void main()
{
	pos = gl_Vertex + vec4(cloud_offset, 1.0);
	gl_Position = ftransform();
	set_fog_coord();
	color = gl_Color;
} 
