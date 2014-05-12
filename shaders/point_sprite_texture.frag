uniform sampler2D tex0;
uniform float min_alpha = 0.0;

void main()
{
	vec4 texel = texture2D(tex0, gl_PointCoord);
	if (texel.a <= min_alpha) discard;
	fg_FragColor = gl_Color*texel;
}
