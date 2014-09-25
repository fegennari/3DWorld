uniform sampler2D color_map, normal_map;
uniform vec4 color_scale     = vec4(1.0);
uniform vec2 normal_tc_off   = vec2(0.0);
uniform vec2 normal_tc_scale = vec2(1.0);

varying vec4 eye_space_pos;
varying vec2 tc;

void main()
{
	vec2 tc_scaled = normal_tc_scale*tc;
	vec4 texel = texture2D(color_map, tc_scaled);
	if (texel.a < 0.75) discard; // transparent
	//if (normal.w == 0.0) discard; // normal not written to (uses nearest filter)
	check_noise_and_maybe_discard(0.0, gl_Color.a);

	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = normalize(fg_NormalMatrix * (2.0*texture2D(normal_map, (tc_scaled + normal_tc_off)).xyz - vec3(1.0)));
	if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse (could use faceforward())
	
	vec3 color = vec3(0.0);
	const bool shadowed = false;
	if (enable_light0) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 0).rgb;
	if (enable_light1) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 1).rgb;
	if (enable_light2) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 2).rgb * calc_light_atten(eye_space_pos, 2);
	fg_FragColor = apply_fog_epos(vec4(clamp(color*color_scale.rgb, 0.0, 1.0)*texel.rgb, 1.0), eye_space_pos);
}
