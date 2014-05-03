varying vec3 dlpos, dl_normal; // world space
varying vec3 normal;
varying vec4 epos; // needed for detail normal map

void main()
{
	setup_texgen_st();
	dl_normal = normalize(fg_Normal);
	vec3 n = fg_NormalMatrix * fg_Normal;
	normal = (no_normalize ? n : normalize(n));
	dlpos  = fg_Vertex.xyz;
	epos   = fg_ModelViewMatrix * fg_Vertex; // needed for detail normal map
	gl_Position   = fg_ftransform();
	gl_FrontColor = fg_Color;
}
