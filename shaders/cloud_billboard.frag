uniform sampler2D tex0;

varying vec2 tc;

void main()
{
	vec2 lum_inv_alpha = texture2D(tex0, tc).rb;
	fg_FragColor = vec4((lum_inv_alpha[0] * gl_Color.rgb), ((1.0 - lum_inv_alpha[1]) * gl_Color.a));
	//fg_FragColor = apply_fog(fg_FragColor);
}
