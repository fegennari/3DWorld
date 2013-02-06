uniform mat4 world_space_mvm;
varying vec4 epos;
varying vec3 normal, world_normal, world_space_pos, vertex;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0;
	world_normal    = gl_Normal;
	normal          = normalize(gl_NormalMatrix * gl_Normal);
	vertex          = gl_Vertex.xyz;
	epos            = gl_ModelViewMatrix * gl_Vertex;
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
	gl_Position     = ftransform();
	gl_FrontColor   = gl_Color;
}

