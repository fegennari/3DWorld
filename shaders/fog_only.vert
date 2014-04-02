void main()
{
	gl_Position = fg_ftransform();
	set_fog_coord(fg_Vertex);
} 
