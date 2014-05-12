uniform sampler2D branch_tex;
uniform float min_alpha = 0.0;
uniform float min_noise = 0.0;
uniform float max_noise = 1.0;

varying float world_space_zval;
varying vec2 tc;

void main()
{
	vec4 texel = texture2D(branch_tex, tc);
	if (texel.a <= min_alpha) discard;
	check_noise_and_maybe_discard(min_noise, max_noise);
	fg_FragColor = apply_fog_scaled(vec4(texel.rgb * gl_Color.rgb, 1.0), world_space_zval);
}
