uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float opacity   = 1.0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	gl_FragColor = apply_fog(gl_Color*texel);
}
