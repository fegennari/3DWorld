uniform sampler2D tex0, alpha_mask_tex;
uniform float min_alpha  = 0.0;
uniform float lum_scale  = 0.0;
uniform float lum_offset = 0.0;
uniform vec4  emission   = vec4(0,0,0,1);

// thse come from bump_map.part
//in vec2 tc;
//in vec4 epos;
//in vec3 eye_norm;

subroutine void postproc_color(); // signature
subroutine(postproc_color) void no_op() {}
subroutine(postproc_color) void apply_burn_mask() {
	if (burn_offset > -1.0) {fg_FragColor = apply_burn_mask(fg_FragColor, tc);} // slow
}
subroutine uniform postproc_color postproc_color_op;

subroutine void maybe_bump_map(inout vec3 normal, inout vec3 light_dir, inout vec3 eye_pos); // signature
subroutine(maybe_bump_map) void  no_bump_map(inout vec3 normal, inout vec3 light_dir, inout vec3 eye_pos) {}
subroutine(maybe_bump_map) void yes_bump_map(inout vec3 normal, inout vec3 light_dir, inout vec3 eye_pos) {
	maybe_apply_bump_map_self_shadowed(normal, light_dir, eye_pos, 1.0);
}
subroutine uniform maybe_bump_map maybe_bump_map_op;

subroutine void do_lighting(inout vec3 color, in vec4 texel, in vec3 n); // signature

vec4 add_ship_light(in vec3 normal, in int i) {
	vec4 epos2     = epos; // copy so that we can modify it
	vec3 light_dir = normalize(fg_LightSource[i].position.xyz - epos2.xyz);
	maybe_bump_map_op(normal, light_dir, epos2.xyz);
	return get_ads_lighting(normal, light_dir, epos2, 1.0, 1.0, gl_Color, fg_LightSource[i]) * calc_light_atten(epos2, i);
}

subroutine(do_lighting) void shadow_only(inout vec3 color, in vec4 texel, in vec3 n) {
	color += add_ship_light(n, 0).rgb; // light 0 is the system light
}
subroutine(do_lighting) void normal_lighting(inout vec3 color, in vec4 texel, in vec3 n) {
	for (int i = 0; i < 8; ++i) {color += add_ship_light(n, i).rgb;}
	// add other emissive term based on texture luminance
	color = mix(color, vec3(1.0), clamp(lum_scale*(texel.r + texel.g + texel.b + lum_offset), 0.0, 1.0));
}
subroutine uniform do_lighting do_lighting_op;

void main()
{
#ifdef ALPHA_MASK_TEX
	if (min_alpha != 0.0) { // this test may or may not help performance
		float alpha_mask = texture(alpha_mask_tex, tc).r;
		if (alpha_mask < min_alpha) discard; // slow
	}
#endif
	vec4 texel = texture(tex0, tc);
	//if (texel.a <= min_alpha) discard; // slow
	vec3 n = (gl_FrontFacing ? normalize(eye_norm) : -normalize(eye_norm)); // two-sided lighting
	vec3 color = emission.rgb;
	do_lighting_op(color, texel, n);
	fg_FragColor = vec4(texel.rgb * clamp(color, 0.0, 1.0), texel.a * gl_Color.a); // use gl_Color alpha directly
	postproc_color_op();
}
