uniform sampler2D tex0, tex1;

void main()
{
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);
	gl_FragColor = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), (texel0.a * texel1.a * gl_Color.a));
}
