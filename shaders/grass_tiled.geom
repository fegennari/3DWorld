varying out vec2 tc;

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same

	for (int i = 0; i < 3; ++i) {
		gl_Position = gl_PositionIn[i];
		tc = gl_TexCoordIn[i][0];
		EmitVertex();
	}
}