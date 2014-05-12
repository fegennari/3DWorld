uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec2 tc;

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;
	fg_FragColor = apply_fog(gl_Color*texel);
}
