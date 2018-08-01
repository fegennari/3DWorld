out vec3 vpos, normal; // world space
out vec3 eye_norm;
out vec4 epos; // declared in bump_map.part.frag, but unused in this shader
out vec2 tc; // declared in bump_map.part.frag, but unused in this shader

void main()
{
	gl_Position = fg_ftransform();
	fg_Color_vf = fg_Color;
	normal   = normalize(fg_Normal);
	eye_norm = normalize(fg_NormalMatrix * normal);
	vpos     = fg_Vertex.xyz;
} 
