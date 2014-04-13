varying vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	gl_Position   = fg_ftransform();
	vec3 normal   = normalize(fg_NormalMatrix * fg_Normal); // eye space
	gl_FrontColor = add_light_comp(normal, 4); // only light 4
}
