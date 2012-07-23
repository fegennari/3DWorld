uniform sampler2D branch_tex, snow_tex, noise_tex;
uniform float min_alpha = 0.0;
uniform float noise_tex_size = 128.0;
uniform float min_noise = 0.0;
uniform float max_noise = 1.0;
uniform float camera_facing_scale = 0.0;

void main()
{
	vec4 texel = texture2D(branch_tex, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	float noise_val = texture2D(noise_tex, gl_FragCoord.xy/noise_tex_size).r;
	if (noise_val < min_noise || noise_val > max_noise) discard;
	vec4 color = texel * vec4(gl_Color.rgb, 1.0);

	if (add_snow) {
		float sw  = (gl_FrontFacing ? 1.0 : camera_facing_scale)*(1.0 - gl_Color.a);
		color.rgb = (1.0 - sw)*color.rgb + sw*0.75*texture2D(snow_tex, gl_TexCoord[0].st).rgb;
	}
	gl_FragColor = apply_fog(color);
}
