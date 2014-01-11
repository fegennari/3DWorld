uniform sampler2D branch_tex;
uniform float min_alpha = 0.0;
uniform float min_noise = 0.0;
uniform float max_noise = 1.0;

varying float world_space_zval;

void main()
{
	vec4 texel = texture2D(branch_tex, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	check_noise_and_maybe_discard(min_noise, max_noise);
	gl_FragColor = apply_fog_scaled(vec4(texel.rgb * gl_Color.rgb, 1.0), world_space_zval);
}
