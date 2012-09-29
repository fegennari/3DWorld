uniform float min_alpha = 0.0;

// clipped eye position, clipped vertex position
varying vec3 eye, vpos, normal; // world space
varying vec4 epos;
varying vec3 eye_norm;

void main()
{
	vec3 norm_normal = normalize(normal);
	vec4 texel;

	if (use_noise_tex) {
		texel = get_texture_val(norm_normal, vpos);
	}
	else {
		texel = lookup_triplanar_texture(vpos, norm_normal, tex0, tex0, tex0) * color0;
	}
	float alpha = gl_Color.a;
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
	add_indir_lighting(lit_color);
	
	// directional light sources with no attenuation (Note: could add other lights later)
	if (enable_light0)  lit_color += add_light_comp_pos_smap(normalize(eye_norm), epos, 0).rgb;
	if (enable_light1)  lit_color += add_light_comp_pos_smap(normalize(eye_norm), epos, 1).rgb;
	if (enable_dlights) lit_color += add_dlights(vpos, norm_normal, eye, gl_FrontMaterial.diffuse.rgb); // dynamic lighting
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif
#ifndef NO_FOG
	color = apply_fog(color); // apply standard fog
#endif
	gl_FragColor = color;
}
