uniform sampler2D tex0, tex1;
uniform float min_alpha = 0.0;

void main()
{
	vec4 texel0  = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1  = texture2D(tex1, gl_TexCoord[1].st);
	if (texel0.a*texel1.a*gl_Color.a < min_alpha) discard;
	gl_FragColor = (texel0 * texel1 * gl_Color); // add fog?
}
