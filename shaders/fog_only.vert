void main()
{
	gl_Position = ftransform();
	set_fog_coord(gl_Vertex);
} 
