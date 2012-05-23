attribute vec3 tangent;

void main()
{
	gl_Position  = gl_Vertex;
	gl_TexCoord[0].xyz = gl_Normal;
	gl_TexCoord[1].xyz = tangent;
	calc_leaf_lighting();
} 
