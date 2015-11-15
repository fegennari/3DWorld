uniform sampler2D tex0;
uniform float opacity = 1.0;

in vec4 epos;
in vec3 normal; // eye space
in vec2 tc;

void main() {
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	vec4 texel   = texture(tex0, tc);
	vec3 normal2 = normalize(normal); // renormalize
	vec3 color   = do_shadowed_lighting(vec4(0.0), epos, normal2, gl_Color, 1.0, 1.0);
	fg_FragColor = apply_fog_epos(texel*vec4(color, 1.0), epos);
}
