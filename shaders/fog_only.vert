void main()
{
	gl_Position = ftransform();
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz);
} 
