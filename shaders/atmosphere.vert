uniform mat4 fg_ViewMatrix;
varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	normal          = normalize(fg_NormalMatrix * fg_Normal);
	epos            = fg_ModelViewMatrix * fg_Vertex;
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position     = fg_ProjectionMatrix * epos;
	gl_FrontColor   = vec4(1.0); // white
}
