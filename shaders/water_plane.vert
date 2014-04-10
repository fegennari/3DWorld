uniform float normal_z = 1.0;

varying vec3 normal;
varying vec4 epos, proj_pos;
varying vec2 tc2;

void main()
{
	setup_texgen_st();
	tc2    = fg_TexCoord;
	normal = gl_NormalMatrix * vec3(0.0, 0.0, normal_z); // not normalized
	epos   = gl_ModelViewMatrix * fg_Vertex;
	proj_pos = fg_ftransform();
	gl_Position = proj_pos;
	gl_FrontColor = fg_Color;
}
