out vec4 epos, proj_pos;
out vec2 tc2;

void main()
{
	setup_texgen_st();
	tc2         = fg_TexCoord;
	epos        = fg_ModelViewMatrix * fg_Vertex;
	proj_pos    = fg_ProjectionMatrix * epos;
	gl_Position = proj_pos;
	fg_Color_vf = fg_Color;
}
