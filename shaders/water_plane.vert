uniform float normal_z = 1.0;

out vec3 normal;
out vec4 epos, proj_pos;
out vec2 tc2;

void main()
{
	setup_texgen_st();
	tc2    = fg_TexCoord;
	normal = fg_NormalMatrix * vec3(0.0, 0.0, normal_z); // not normalized
	epos   = fg_ModelViewMatrix * fg_Vertex;
	proj_pos    = fg_ProjectionMatrix * epos;
	gl_Position = proj_pos;
	fg_Color_vf = fg_Color;
}
