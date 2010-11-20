void set_fog()
{
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz);
}
