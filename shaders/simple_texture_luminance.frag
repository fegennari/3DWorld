uniform sampler2D tex0;

in vec2 tc;

void main()
{
	fg_FragColor = gl_Color*vec4(texture(tex0, tc).rrr, 1.0);
}
