attribute vec3 tangent;

void main()
{
	gl_Position = gl_Vertex;
	gl_TexCoord[6].xyz = normalize(gl_Normal);
	gl_TexCoord[7].xyz = tangent;
	calc_leaf_lighting();
} 
