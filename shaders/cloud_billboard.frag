uniform sampler2D tex0;

void main()
{
	vec2 lum_inv_alpha = texture2D(tex0, gl_TexCoord[0].st).rb;
	vec4 color = vec4((lum_inv_alpha[0] * gl_Color.rgb), ((1.0 - lum_inv_alpha[1]) * gl_Color.a));
	//gl_FragColor = apply_fog(color);
	gl_FragColor = color;
}
