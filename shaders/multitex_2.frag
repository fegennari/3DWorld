uniform sampler2D tex0, tex1;
uniform float min_alpha = 0.0;

in vec2 tc;

void main()
{
	vec4 texel0  = texture(tex0, tc);
	vec4 texel1  = texture(tex1, tc);
	if (texel0.a*texel1.a*gl_Color.a < min_alpha) discard;
	fg_FragColor = (texel0 * texel1 * gl_Color); // add fog?
}
