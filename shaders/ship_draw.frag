uniform sampler2D tex0, alpha_mask_tex;
uniform float min_alpha  = 0.0;
uniform float lum_scale  = 0.0;
uniform float lum_offset = 0.0;
uniform vec4  emission   = vec4(0,0,0,1);

varying vec2 tc;
varying vec4 epos;
varying vec3 normal; // eye space

void main()
{
#ifdef ALPHA_MASK_TEX
	if (min_alpha != 0.0) { // this test may or may not help performance
		float alpha_mask = texture2D(alpha_mask_tex, tc).r;
		if (alpha_mask < min_alpha) discard; // slow
	}
#endif
	vec4 texel = texture2D(tex0, tc);
	//if (texel.a <= min_alpha) discard; // slow
	vec3 n = (gl_FrontFacing ? normalize(normal) : -normalize(normal)); // two-sided lighting
	vec3 color = emission.rgb;
#ifdef SHADOW_ONLY_MODE
	color += add_pt_light_comp(n, epos, 0).rgb; // light 0 is the system light
#else
	for (int i = 0; i < 8; ++i) {color += add_pt_light_comp(n, epos, i).rgb;}
	// add other emissive term based on texture luminance
	color = mix(color, vec3(1.0), clamp(lum_scale*(texel.r + texel.g + texel.b + lum_offset), 0.0, 1.0));
#endif
	fg_FragColor = vec4(texel.rgb * clamp(color, 0.0, 1.0), texel.a * gl_Color.a); // use gl_Color alpha directly
#ifdef BURN_MASK_TEX
	if (burn_offset > -1.0) {fg_FragColor = apply_burn_mask(fg_FragColor, tc);} // slow, conditional may or may not help
#endif
}
