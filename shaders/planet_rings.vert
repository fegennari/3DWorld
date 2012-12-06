uniform mat4 world_space_mvm;
varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0; // remove later
	normal          = normalize(gl_NormalMatrix * gl_Normal);
	epos            = gl_ModelViewMatrix * gl_Vertex;
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
	gl_Position     = ftransform();
	gl_FogFragCoord = length(epos.xyz);
	gl_FrontColor   = gl_Color;
}
