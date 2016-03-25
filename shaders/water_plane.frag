#include <fresnel.part>

uniform float normal_z = 1.0;
uniform sampler2D reflection_tex, water_normal_tex, height_tex, noise_tex, deep_water_normal_tex, foam_tex, shadow_tex;
uniform vec4 water_color, reflect_color;
uniform float noise_time, wave_time, wave_amplitude, water_plane_z, water_green_comp, reflect_scale, mesh_z_scale, rain_intensity;

in vec4 epos, proj_pos;
in vec2 tc, tc2;
#ifdef TESS_MODE
in float water_zval, wave_dz;
in vec3 normal;
#endif

vec3 water_normal_lookup(in vec2 wtc) {
	return 2.0*(texture(water_normal_tex, 0.5*wtc).rgb - 0.5);
}
vec3 get_wave_normal(in vec2 wtc) {
	float ntime = 2.0*abs(fract(0.005*wave_time) - 0.5);
	return mix(water_normal_lookup(wtc), water_normal_lookup(wtc + vec2(0.3, 0.6)), ntime);
}

vec3 deep_water_normal_lookup(in vec2 wtc) {
	return 2.0*(texture(deep_water_normal_tex, 1.04*wtc).rgb - 0.5);
}
vec4 get_norm_foam_val(vec3 n) {
	return vec4(n, clamp(3.0*(0.5 - dot(n, vec3(0,0,1))), 0.0, 1.0)); // use dot product with +z
}
vec4 get_deep_wave_normal(in vec2 wtc) {
	float ntime = 2.0*abs(fract(0.011111*wave_time) - 0.5);
	vec3 n1 = deep_water_normal_lookup(wtc);
	vec3 n2 = deep_water_normal_lookup(wtc + vec2(0.5, 0.5));
	// reduce the pulsing effect by increasing normal map amplitude when two normal maps are averaged together
	float mag = sqrt(2.0 - 2.0*abs(ntime - 0.5)); 
	return mix(get_norm_foam_val(n1), get_norm_foam_val(n2), ntime)*mag;
}


void main() {
#ifdef USE_WATER_DEPTH // else we assume water is neither too shallow nor too deep
	float mesh_z = texture(height_tex, tc2).r;
#ifdef TESS_MODE
	float depth  = water_zval - mesh_z;
#else
	float depth  = water_plane_z - mesh_z;
#endif
	if (depth <= 0.0) discard;
#endif // USE_WATER_DEPTH

#ifdef WRITE_DEPTH_ONLY
	fg_FragColor = vec4(normal, 1.0); // Note: for some reason, we can have depth artifacts if we don't use the normal
#else // !WRITE_DEPTH_ONLY

#ifdef TESS_MODE
	vec3 norm = normalize(normal);
#else
	vec3 norm = normal_z*fg_NormalMatrix[2]; // eye space (+/- z in world space)
#endif // TESS_MODE
	vec2 ripple    = vec2(0,0);
	vec3 add_color = vec3(0);
	float foam_amt = 0.0;

	if (add_noise) { // for rain
		vec3 wave_n = get_wave_normal(fract(4.61*noise_time)*3.0*proj_pos.xy/proj_pos.w);
		norm        = normalize(norm + 0.06*fg_NormalMatrix*wave_n);
		ripple     += 0.025*wave_n.xy;
		vec2 wtc    = 20.0*tc2 + fract(vec2(1.7*noise_time, 1.6*noise_time));
		add_color  += vec3(1.0) * (1.0 - texture(noise_tex, wtc).r); // Note that the texture is white with blue dots
	}
	vec3 light_norm = norm;
	float green_scale = 0.0;

	if (add_waves) { // deep water waves are the default
		// calculate ripple adjustment of normal and reflection based on scaled water normal map texture
		vec4 norm_fa = get_deep_wave_normal(tc);
		vec3 wave_n  = 1.25*wave_amplitude*normalize(norm_fa.xyz);
#ifdef USE_WATER_DEPTH
		// deep water waves shouldn't move (much) with the wind, but that would require another set of TCs, texgen matrix, etc.
		float deep_wave_scale = clamp((0.8*depth*mesh_z_scale - 0.2), 0.0, 1.0);
		wave_n = mix(wave_amplitude*get_wave_normal(tc), wave_n, deep_wave_scale);
		if (use_foam) {foam_amt = deep_wave_scale*norm_fa.w;}
#endif // USE_WATER_DEPTH
		vec3 wave_n_eye = fg_NormalMatrix * wave_n;
		if (reflections) {green_scale += 1.0*(1.0 - abs(dot(norm, normalize(wave_n_eye))));} // add green to sides of waves
		light_norm = normalize(norm + wave_n_eye);
		norm       = normalize(norm + 0.1*wave_n_eye); // lower scale for fresnel term
		ripple    += 0.05*wave_n.xy;
	}
	if (add_rain) {
		float ripple_mag = clamp(2.0*(1.0 - 0.25*length(epos.xyz)), 0.0, 1.0);

		if (ripple_mag > 0.0) {
			vec3 ripple_norm = get_ripple_normal(12.0*tc, 0.2*wave_time, rain_intensity);
			ripple_norm.z = 0.0; // clear z-component so that the non-rippled parts of the texture don't affect the wave normals
			ripple_norm = 2.0*ripple_mag * (fg_NormalMatrix * ripple_norm); // in eye space, not normalized
			light_norm  = normalize(ripple_norm + light_norm);
			norm        = normalize(ripple_norm + norm);
		}
	}

	// calculate lighting
	vec3 epos_n          = normalize(epos.xyz);
	float cos_view_angle = abs(dot(epos_n, norm));
	vec4 color           = water_color;
#ifdef USE_WATER_DEPTH
	if (!is_lava) {color.a *= mix(1.0, clamp(20.0*depth, 0.0, 1.0), min(1.0, 2.5*cos_view_angle));} // blend to alpha=0 near the shore
#endif

#ifdef ENABLE_WATER_SHADOWS // looks okay for shallow water, but bad for deep water since we get underwater shadows on the water surface
	vec3 shadow  = texture(shadow_tex, tc2).rgb; // {mesh_shadow, tree_shadow, ambient_occlusion}
	float ascale = shadow.b;
	float dscale = min(shadow.r, shadow.g);
#else
	float ascale = 1.0;
	float dscale = 1.0;
#endif
	vec4 lighting = vec4(0,0,0,1);
	if (is_lava) {lighting += vec4(0.3, 0.05, 0.0, 0.0);} // add emissive light

	// Note: when drawing the mesh we use two passes, one with shadow maps enabled and one without
	// for water we only have one pass, and the shadow map texture may not be bound/valid in all cases,
	// but we know that it's valid when smap_scale > 0.0, so we can/must use this test to control the shadow map texture lookup
	lighting.rgb += do_shadowed_lighting(vec4(0.0), epos, light_norm, gl_Color, ascale, dscale); // Note: vertex parameter is unused
	
#ifdef TESS_MODE
	float subsurf_scatter = max(0.55*(1.0 - cos_view_angle)*wave_dz, 0.0);
	green_scale += 0.2*pow(subsurf_scatter, 4.0);
#endif
	// add some green at shallow view angles
	green_scale += (1.0 - cos_view_angle);
	color = mix(color, vec4(0.0, 1.0, 0.5, color.a), min(1.0, water_green_comp*green_scale));

	if (reflections) { // calculate reflections
		float reflect_w  = reflect_scale*get_fresnel_reflection(-epos_n, norm, 1.0, 1.333);
		//float reflect_w  = reflect_scale*get_approx_fresnel_reflection(-epos_n, norm, 1.0, 1.333);
		vec2 ref_tex_st  = clamp(0.5*proj_pos.xy/proj_pos.w + 0.3*ripple + vec2(0.5, 0.5), 0.0, 1.0);
		vec4 reflect_tex = vec4(texture(reflection_tex, ref_tex_st).rgb, 1.0);
		color = mix(color, reflect_color * reflect_tex, reflect_w);
	}
#ifdef USE_WATER_DEPTH
	if (use_foam) {
		// foam texture near shore (what about is_lava?)
		foam_amt += 0.5*clamp(10.0*depth-0.2, 0.0, 1.0)*(1.0 - clamp(5.0*depth, 0.0, 1.0));
		color = mix(color, texture(foam_tex, 25.0*tc), foam_amt);
	}
#endif
	// determine final color with fog
	color.rgb   += add_color;
	fg_FragColor = vec4(color.rgb * lighting.rgb, color.a * gl_Color.a); // use gl_Color alpha directly
	fg_FragColor = apply_fog_epos(fg_FragColor, epos);
#endif // WRITE_DEPTH_ONLY
}
