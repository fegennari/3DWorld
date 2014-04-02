uniform mat4 world_space_mvm;
varying vec3 dlpos, dl_normal; // world space
varying vec3 normal;

void main()
{
	setup_texgen0();
	dl_normal = normalize(fg_Normal);
	vec3 n = gl_NormalMatrix * fg_Normal;
	normal = (no_normalize ? n : normalize(n));
	dlpos  = fg_Vertex.xyz;
	gl_Position = fg_ftransform();
	gl_FrontColor = fg_Color;
}
