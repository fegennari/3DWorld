uniform sampler2D tex0, tex1;
uniform float min_alpha = 0.0;

varying vec2 tc;

void main()
{
	vec4 texel0  = texture2D(tex0, tc);
	vec4 texel1  = texture2D(tex1, tc);
	if (texel0.a*texel1.a*gl_Color.a < min_alpha) discard;
	gl_FragColor = (texel0 * texel1 * gl_Color); // add fog?
}
