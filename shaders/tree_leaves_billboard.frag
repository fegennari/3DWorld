uniform sampler2D tex0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a == 0.0 || (texel.r + texel.g + texel.b) < 0.5) discard; // background is black
	gl_FragColor = apply_fog(gl_Color*texel);
}
