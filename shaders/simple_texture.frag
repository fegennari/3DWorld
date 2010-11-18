uniform sampler2D tex;
uniform float min_alpha = 0.0;

void main()
{
	vec4 texel = texture2D(tex, gl_TexCoord[0].st);
	if (texel.a < min_alpha) discard;
	gl_FragColor = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
}
