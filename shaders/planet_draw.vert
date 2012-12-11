varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0;
	normal          = normalize(gl_NormalMatrix * gl_Normal);
	world_space_pos = gl_Vertex.xyz;
	epos            = gl_ModelViewMatrix * gl_Vertex;
	gl_Position     = ftransform();
	gl_FogFragCoord = length(epos.xyz);
	gl_FrontColor   = gl_Color;
}

