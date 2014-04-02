uniform mat4 world_space_mvm;
varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;

void main()
{
	tc              = fg_TexCoord;
	normal          = normalize(gl_NormalMatrix * fg_Normal);
	vertex          = fg_Vertex.xyz;
	vec4 epos       = gl_ModelViewMatrix * fg_Vertex;
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
	gl_Position     = fg_ftransform();
	gl_FrontColor   = vec4(1.0); // always white - color will come from the texture
}

