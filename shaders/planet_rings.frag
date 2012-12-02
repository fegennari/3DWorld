uniform sampler2D tex0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	gl_FragColor = apply_fog(gl_Color * texel); // FIXME: make more interesting
}
