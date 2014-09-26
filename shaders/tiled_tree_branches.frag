uniform sampler2D tex0;
uniform float opacity = 1.0;

in vec4 epos;
in vec3 normal;
in vec2 tc;

void main()
{
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	vec4 texel   = texture2D(tex0, tc);
	vec3 normal2 = normalize(normal); // renormalize
	vec3 color   = vec3(0.0);
	if (enable_light0) color += add_light_comp_pos0(normal2, epos).rgb; // sun
	if (enable_light1) color += add_light_comp_pos1(normal2, epos).rgb; // moon
	if (enable_light2) color += add_light_comp_pos (normal2, epos, 2).rgb * calc_light_atten(epos, 2); // lightning
	fg_FragColor = apply_fog_epos(texel*vec4(color, 1.0), epos);
}
