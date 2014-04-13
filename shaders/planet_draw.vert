uniform mat4 fg_ViewMatrix;
varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;

void main()
{
	tc              = fg_TexCoord;
	normal          = normalize(fg_NormalMatrix * fg_Normal);
	vertex          = fg_Vertex.xyz;
	vec4 epos       = fg_ModelViewMatrix * fg_Vertex;
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position     = fg_ProjectionMatrix * epos;
	gl_FrontColor   = vec4(1.0); // always white - color will come from the texture
}

