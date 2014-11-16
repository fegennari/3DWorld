uniform sampler2D tex0;
uniform float min_alpha = 0.0;

in vec2 tc;

void main()
{
	if (texture(tex0, tc).r <= min_alpha) discard;
	fg_FragColor = gl_Color;
}
