uniform sampler2D tex;

void main()
{
	vec4 texel = texture2D(tex,gl_TexCoord[0].st);
	gl_FragColor = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
}
