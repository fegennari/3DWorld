uniform float min_alpha = 0.0;
uniform vec4 color0 = vec4(1,1,1,1);
uniform vec4 color1 = vec4(1,1,1,1);
uniform sampler2D tex0, tex1;

// clipped eye position, clipped vertex position
varying vec3 eye, vpos, normal; // world space
varying vec4 epos;
varying vec3 eye_norm;

void main()
{
	vec3 norm_normal = normalize(normal);
	vec4 texel0 = lookup_triplanar_texture(vpos, norm_normal, tex0, tex0, tex0) * color0;
	vec4 texel;

	if (use_noise_tex) {
		vec4 texel1 = lookup_triplanar_texture(vpos, norm_normal, tex1, tex1, tex1) * color1;
		// interpolate between the two texture/color pairs using a random noise function
		texel = mix(texel0, texel1, procedural_eval(vpos));
	}
	else {
		texel = texel0;
	}
	float alpha = gl_Color.a;
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
	add_indir_lighting(lit_color);
	
	// directional light sources with no attenuation (Note: could add other lights later)
	if (enable_light0)  lit_color += add_light_comp_pos_smap(normalize(eye_norm), epos, 0).rgb;
	if (enable_light1)  lit_color += add_light_comp_pos_smap(normalize(eye_norm), epos, 1).rgb;
	if (enable_dlights) lit_color += gl_FrontMaterial.diffuse.rgb * add_dlights(vpos, norm_normal, eye); // dynamic lighting
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif
#ifndef NO_FOG
	color = apply_fog(color); // apply standard fog
#endif
	gl_FragColor = color;
}
