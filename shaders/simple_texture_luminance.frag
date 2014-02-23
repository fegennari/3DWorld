uniform sampler2D tex0;

varying vec2 tc;

void main()
{
	gl_FragColor = gl_Color*vec4(texture2D(tex0, tc).rrr, 1.0);
}
