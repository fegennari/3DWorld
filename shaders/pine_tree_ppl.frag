uniform sampler2D branch_tex;
uniform float min_alpha = 0.0;
uniform float min_noise = 0.0;
uniform float max_noise = 1.0;

in float world_space_zval;
in vec2 tc;
in vec4 epos;
in vec3 normal; // eye space

void main() {
	vec4 texel = texture(branch_tex, tc);
	if (texel.a <= min_alpha) discard;
#ifndef NO_NOISE
	check_noise_and_maybe_discard(min_noise, max_noise);
#endif
	vec3 color   = do_shadowed_lighting(vec4(0.0), epos, normal, gl_Color, 1.0, 1.0);
	fg_FragColor = apply_fog_scaled(vec4(texel.rgb * color, 1.0), world_space_zval);
}
