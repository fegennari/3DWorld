uniform sampler2D tex0;
uniform float min_alpha = 0.0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	gl_FragColor = apply_fog(gl_Color*texel);
}
