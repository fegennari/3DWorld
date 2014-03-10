varying vec3 normal;
varying vec4 epos, proj_pos;
varying vec2 tc, tc2;

uniform sampler2D reflection_tex, water_normal_tex, height_tex, noise_tex, deep_water_normal_tex;
uniform vec4 water_color, reflect_color;
uniform float noise_time, wave_time, wave_amplitude, water_plane_z, water_green_comp, reflect_scale, mesh_z_scale;

vec3 water_normal_lookup(in vec2 wtc) {
	return 2.0*(texture2D(water_normal_tex, 0.5*wtc).rgb - 0.5);
}
vec3 get_wave_normal(in vec2 wtc) {
	float ntime = 2.0*abs(fract(0.005*wave_time) - 0.5);
	return mix(water_normal_lookup(wtc), water_normal_lookup(wtc + vec2(0.3, 0.6)), ntime);
}

vec3 deep_water_normal_lookup(in vec2 wtc) {
	return 2.0*(texture2D(deep_water_normal_tex, 1.04*wtc).rgb - 0.5);
}
vec4 get_norm_foam_val(vec3 n) {
	return vec4(n, clamp(3.0*(0.5 - dot(n, vec3(0,0,1))), 0.0, 1.0)); // use dot product with +z
}
vec4 get_deep_wave_normal(in vec2 wtc) {
	float ntime = 2.0*abs(fract(0.011*wave_time) - 0.5);
	vec3 n1 = deep_water_normal_lookup(wtc);
	vec3 n2 = deep_water_normal_lookup(wtc + vec2(0.5, 0.5));
	return mix(get_norm_foam_val(n1), get_norm_foam_val(n2), ntime);
}

void main()
{
#ifdef USE_WATER_DEPTH // else we assume water is neither too shallow nor too deep
	float mesh_z = texture2D(height_tex, tc2).r;
	float depth  = water_plane_z - mesh_z;
	if (depth <= 0.0) discard;
#endif
	vec3 norm   = normalize(normal); // renormalize
	vec2 ripple = vec2(0,0);
	vec3 add_color = vec3(0);
	float foam_amt = 0.0;

	if (add_noise) { // for rain
		vec3 wave_n = get_wave_normal(fract(4.61*noise_time)*3.0*proj_pos.xy/proj_pos.w);
		norm        = normalize(norm + 0.06*gl_NormalMatrix*wave_n);
		ripple     += 0.025*wave_n.xy;
		vec2 wtc    = 20.0*tc2 + vec2(2.7*noise_time, 2.6*noise_time);
		add_color  += vec3(1) * (1.0 - texture2D(noise_tex, wtc).r); // Note that the texture is white with blue dots
	}
	vec3 light_norm = norm;
	float green_scale = 0.0;

	if (add_waves) {
		// calculate ripple adjustment of normal and reflection based on scaled water normal map texture
		vec3 wave_n = wave_amplitude*get_wave_normal(tc);
#ifdef USE_WATER_DEPTH
		if (deep_water_waves) {
			// deep water waves shouldn't move (much) with the wind, but that would require another set of TCs, texgen matrix, etc.
			vec4 norm_fa = get_deep_wave_normal(tc);
			float deep_wave_scale = clamp((0.8*depth*mesh_z_scale - 0.2), 0.0, 1.0);
			wave_n   = mix(wave_n, 1.25*wave_amplitude*normalize(norm_fa.xyz), deep_wave_scale);
			foam_amt = deep_wave_scale*norm_fa.w;
		}
#endif
		vec3 wave_n_eye = gl_NormalMatrix * wave_n;
		if (reflections) {green_scale += 0.8*(1.0 - abs(dot(norm, normalize(wave_n_eye))));} // add green to sides of waves (though this increases shader time)
		light_norm = normalize(norm + wave_n_eye);
		norm       = normalize(norm + 0.1*wave_n_eye); // lower scale for fresnel term
		ripple    += 0.05*wave_n.xy;
	}

	// calculate lighting
	vec3 epos_n   = normalize(epos.xyz);
	float cos_view_angle = abs(dot(epos_n, norm));
	vec4 color    = water_color;
#ifdef USE_WATER_DEPTH
	color.a      *= mix(1.0, clamp(20.0*depth, 0.0, 1.0), min(1.0, 2.5*cos_view_angle)); // blend to alpha=0 near the shore
#endif
	vec4 lighting = vec4(0,0,0,1);
	if (enable_light0) {lighting += add_light_comp_pos0(light_norm, epos);}
	if (enable_light1) {lighting += add_light_comp_pos1(light_norm, epos);}
	
	// add some green at shallow view angles
	green_scale += (1.0 - cos_view_angle);
	color = mix(color, vec4(0.0, 1.0, 0.5, color.a), water_green_comp*min(1.0, green_scale));
	color = mix(color, vec4(1.0), foam_amt);

	if (reflections) { // calculate reflections
		float reflect_w  = reflect_scale*get_fresnel_reflection(-epos_n, norm, 1.0, 1.333);
		vec2 ref_tex_st  = clamp(0.5*proj_pos.xy/proj_pos.w + 0.3*ripple + vec2(0.5, 0.5), 0.0, 1.0);
		vec4 reflect_tex = vec4(texture2D(reflection_tex, ref_tex_st).rgb, 1.0);
		color = mix(color, reflect_color * reflect_tex, reflect_w);
	}

	// determine final color with fog
	color.rgb += add_color;
	vec4 frag_color = vec4(color.rgb * lighting.rgb, color.a * gl_Color.a); // use gl_Color alpha directly
	gl_FragColor = apply_fog_epos(frag_color, epos);
}
