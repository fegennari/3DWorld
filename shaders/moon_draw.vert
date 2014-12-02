out vec2 tc;
out vec3 normal;

void main()
{
	tc          = fg_TexCoord;
	gl_Position = fg_ftransform();
	normal      = normalize(fg_NormalMatrix * fg_Normal); // eye space
	fg_Color_vf = gl_Color;
}
