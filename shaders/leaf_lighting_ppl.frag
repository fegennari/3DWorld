uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float opacity   = 1.0;
uniform vec4 color_scale = vec4(1.0);

in vec2 tc;
in vec4 epos;
in vec3 normal; // eye space
in vec3 ws_pos;
in vec3 ws_normal;

void main() {
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;
#ifdef ENABLE_OPACITY
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
#endif
	vec3 ws_normal2 = ((dot(normal, epos.xyz) > 0.0) ? -1.0 : 1.0)*normalize(ws_normal);
	vec3 color      = get_tree_leaf_lighting(epos, normal, ws_pos, ws_normal2, 1.0);
	fg_FragColor = texel*vec4(min(2.0*gl_Color.rgb, clamp(color*color_scale.rgb, 0.0, 1.0)), 1.0); // limit lightning color
#ifndef NO_FOG
	fg_FragColor = apply_fog(fg_FragColor);
#endif
}
