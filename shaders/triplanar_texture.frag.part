uniform sampler2D tex0;
uniform float min_alpha = 0.0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	gl_FragColor = apply_fog(color);
}
