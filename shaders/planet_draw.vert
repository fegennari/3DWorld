uniform mat4 world_space_mvm;
varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;

void main()
{
	tc              = gl_MultiTexCoord0;
	normal          = normalize(gl_NormalMatrix * gl_Normal);
	vertex          = gl_Vertex.xyz;
	vec4 epos       = gl_ModelViewMatrix * gl_Vertex;
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
	gl_Position     = ftransform();
	gl_FrontColor   = vec4(1.0); // always white - color will come from the texture
}

