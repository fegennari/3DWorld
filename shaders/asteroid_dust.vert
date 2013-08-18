uniform mat4 world_space_mvm;

varying vec3 world_space_pos;
varying vec4 epos;

void main()
{
	epos          = gl_ModelViewMatrix * gl_Vertex;
#ifdef ENABLE_SHADOWS
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
#endif
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
}
