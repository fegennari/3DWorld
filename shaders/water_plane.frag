varying vec3 normal;
varying vec4 epos, proj_pos;
uniform sampler2D water_tex, ripple_tex, noise_tex, reflection_tex;
uniform vec4 water_color, reflect_color;
uniform float noise_time;

void main()
{
	vec2 st     = gl_TexCoord[0].st;
	vec4 color  = texture2D(water_tex, st) * water_color;
	vec3 norm   = normalize(normal); // renormalize
	vec2 ripple = vec2(0,0);

	if (add_noise) {
		vec2 st2    = 3.21*proj_pos.xy/proj_pos.w;
		ripple     += 0.1*vec2(texture2D(noise_tex, (1.1*st2 + 12.34*noise_time)).g, texture2D(ripple_tex, (st2 + vec2(0.5,0.5) + 43.21*noise_time)).g) - 0.05;
		norm        = normalize(norm + gl_NormalMatrix*vec3(ripple, 0));
	}
	if (add_waves) {
		// calculate ripple adjustment of normal and reflection based on scaled water texture
		vec2 wave_d = vec2(texture2D(ripple_tex, 11.0*st).g, texture2D(ripple_tex, 10.0*st+vec2(0.5,0.5)).g) - 0.575;
		ripple     += wave_d;
		vec2 delta_n= clamp(2.0*wave_d, 0, 1);
		delta_n    += 0.05*(vec2(gen_cloud_alpha_norm(200*st), gen_cloud_alpha_norm(207*st + vec2(0.3, 0.6))) - 0.5);
		norm        = normalize(norm + gl_NormalMatrix*vec3(delta_n, 0));
	}

	// calculate lighting
	vec4 lighting = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) lighting += add_light_comp_pos(norm, epos, 0);
	if (enable_light1) lighting += add_light_comp_pos(norm, epos, 1);

	if (reflections) {
		vec3 epos_n = normalize(epos.xyz);

		// add some green at shallow view angles
		color = mix(color, vec4(0.0, 1.0, 0.5, color.a), 0.2*(1.0 - abs(dot(epos_n, norm))));

		// calculate reflections
		float reflect_w  = get_fresnel_reflection(-epos_n, norm, 1.0, 1.333);
		vec2 ref_tex_st  = clamp(0.5*proj_pos.xy/proj_pos.w + 0.3*ripple + vec2(0.5, 0.5), 0.0, 1.0);
		vec4 reflect_tex = vec4(texture2D(reflection_tex, ref_tex_st).rgb, 1.0);
		color = mix(color, reflect_color * reflect_tex, reflect_w);
	}

	// determine final color with fog
	vec4 frag_color = vec4(color.rgb * lighting.rgb, color.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	gl_FragColor = apply_fog(frag_color);
}
