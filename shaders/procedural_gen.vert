out vec3 vpos, normal; // world space
out vec3 eye_norm;

void main()
{
	gl_Position = fg_ftransform();
	fg_Color_vf = fg_Color;
	normal   = normalize(fg_Normal);
	eye_norm = normalize(fg_NormalMatrix * normal);
	vpos     = fg_Vertex.xyz;
} 
