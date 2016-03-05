#include <fresnel.part>

uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta, step_delta_shadow;
uniform sampler2D tex0;
uniform sampler3D wet_noise_tex;
uniform sampler2D sky_zval_tex;
uniform float min_alpha       = 0.0;
uniform float water_depth     = 0.0;
uniform float emissive_scale  = 0.0;
uniform float smoke_const_add = 0.0;
uniform float smoke_noise_mag = 0.0;
uniform vec3 smoke_color, sphere_center;
uniform vec3 sun_pos; // used for dynamic smoke shadows line clipping
uniform vec3 fog_time;
uniform float light_atten = 0.0, refract_ix = 1.0;
uniform float metalness   = 0.0;
uniform float cube_bb[6], sphere_radius;
uniform float depth_trans_bias, clip_plane_z, ripple_time, rain_intensity, reflectivity, snow_cov_amt;
uniform float reflect_plane_ztop, reflect_plane_zbot, winding_normal_sign;
uniform vec4 emission = vec4(0,0,0,1);

//in vec3 vpos, normal; // world space, come from indir_lighting.part.frag
// epos, eye_norm, and tc come from bump_map.frag
// camera_pos comes from dynamic_lighting.part

// Note: at most one of these should be enabled
#ifdef ENABLE_REFLECTIONS
uniform sampler2D reflection_tex;
#endif
#ifdef ENABLE_CUBE_MAP_REFLECT
uniform samplerCube reflection_tex;
uniform float cube_map_near_clip = 1.0;
uniform vec3 cube_map_center     = vec3(0.0);
#endif

const float SMOKE_SCALE = 0.25;

// Note: dynamic point lights use reflection vector for specular, and specular doesn't move when the eye rotates
//       global directional lights use half vector for specular, which seems to be const per pixel, and specular doesn't move when the eye translates
#define ADD_LIGHT(i) lit_color += add_pt_light_comp(n, epos, i).rgb

vec3 add_light0(in vec3 n, in float normal_sign, in vec4 base_color) {
	vec3 light_dir = normalize(fg_LightSource[0].position.xyz - epos.xyz); // Note: could drop the -epos.xyz for a directional light
	float nscale   = 1.0;
#ifdef USE_SHADOW_MAP
	float n_dot_l = dot(n, light_dir);
	nscale = ((n_dot_l > 0.0) ? min(1.0, 10.0*n_dot_l)*get_shadow_map_weight_light0(epos, n) : 0.0); // back-facing test
#endif

#ifdef DYNAMIC_SMOKE_SHADOWS
	pt_pair res = clip_line(vpos, sun_pos, smoke_bb);
	vec3 lpos0  = res.v1;
	vec3 vposl  = res.v2;

	if (lpos0 != vposl && nscale > 0.0) {
		vec3 dir      = vposl - lpos0;
		vec3 pos      = (lpos0 - scene_llc)/scene_scale;
		float nsteps  = length(dir)/step_delta_shadow;
		int num_steps = 1 + min(100, int(nsteps)); // round up
		float sd_ratio= max(1.0, nsteps/100.0);
		vec3 delta    = normalize(dir)*sd_ratio*step_delta_shadow/scene_scale;
		float step_weight  = fract(nsteps);
		float smoke_sscale = SMOKE_SCALE*step_delta_shadow/half_dxy;
	
		// smoke volume iteration using 3D texture, light0 to vposl
		for (int i = 0; i < num_steps; ++i) {
			float smoke = smoke_sscale*texture(smoke_and_indir_tex, pos.zxy).a*step_weight*sd_ratio;
			nscale     *= (1.0 - smoke);
			pos        += delta*step_weight; // should be in [0.0, 1.0] range
			step_weight = 1.0;
		}
	}
#endif // DYNAMIC_SMOKE_SHADOWS
	return add_light_comp_pos_scaled_light(nscale*n, epos, 1.0, 1.0, base_color, fg_LightSource[0], normal_sign).rgb;
}

vec3 add_light1(in vec3 n, in float normal_sign, in vec4 base_color) {
#ifdef USE_SHADOW_MAP
	if (use_shadow_map) {n *= get_shadow_map_weight_light1(epos, n);}
#endif
	return add_light_comp_pos_scaled_light(n, epos, 1.0, 1.0, base_color, fg_LightSource[1], normal_sign).rgb;
}

void add_smoke_contrib(in vec3 eye_c, in vec3 vpos_c, inout vec4 color) {
	vec3 dir      = eye_c - vpos_c;
	vec3 norm_dir = normalize(dir); // used for dlights
	vec3 pos      = (vpos_c - scene_llc)/scene_scale;
	float nsteps  = length(dir)/step_delta;
	int num_steps = 1 + min(1000, int(nsteps)); // round up
	float sd_ratio= max(1.0, nsteps/1000.0);
	vec3 delta    = sd_ratio*dir/(nsteps*scene_scale);
	float step_weight  = fract(nsteps);
	float smoke_sscale = SMOKE_SCALE*step_delta/half_dxy;
#ifdef SMOKE_SHADOW_MAP
	vec4 cur_epos   = fg_ModelViewMatrix * vec4(vpos_c, 1.0);
	vec3 epos_delta = fg_NormalMatrix * (delta * scene_scale); // eye space pos/delta
#endif

	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec4 tex_val = texture(smoke_and_indir_tex, pos.zxy); // rgba = {color.rgb, smoke}
#ifdef SMOKE_DLIGHTS
		if (enable_dlights) { // dynamic lighting
			vec3 dl_pos  = pos*scene_scale + scene_llc;
			add_dlights_bm_scaled(tex_val.rgb, dl_pos, norm_dir, vec3(1.0), 0.0, 1.0, 0.0); // normal points from vertex to eye, override bump mapping, color is applied later
		}
#endif // SMOKE_DLIGHTS
#ifdef SMOKE_SHADOW_MAP
#ifdef USE_SHADOW_MAP
		if (enable_light0) {
			const float smoke_albedo = 0.9;
			tex_val.rgb  += smoke_albedo * get_shadow_map_weight_light0_no_bias(cur_epos) * fg_LightSource[0].diffuse.rgb;
			cur_epos.rgb += epos_delta;
		}
#endif // USE_SHADOW_MAP
#endif // SMOKE_SHADOW_MAP

#ifdef SMOKE_NOISE
		tex_val.a += smoke_noise_mag * max(0.0, gen_cloud_alpha_time((pos * vec3(1.0, 1.0, 0.7)), fog_time)); // stretch in z
#endif
		float smoke  = smoke_sscale*(tex_val.a + smoke_const_add)*step_weight*sd_ratio;
		float alpha  = (keep_alpha ? color.a : ((color.a == 0.0) ? smoke : 1.0));
		float mval   = ((!keep_alpha && color.a == 0.0) ? 1.0 : smoke);
		color        = mix(color, vec4((tex_val.rgb * smoke_color), alpha), mval);
		pos         += delta*step_weight; // should be in [0.0, 1.0] range
		step_weight  = 1.0;
	} // for i
}

float get_water_snow_coverage() {
	vec3 pos    = (vpos - scene_llc)/scene_scale;
	float cmp   = vpos.z + 0.1*half_dxy;
	float delta = texture(sky_zval_tex, pos.xy).r - cmp; // incorrectly interpolated, but smooth
	return clamp((1.0 - 5.0*delta/half_dxy), 0.0, 1.0);
}

float get_puddle_val(in float wetness) {

	float wet_val = 0.0;
	float freq    = 1.0;

	for (int i = 0; i < 4; ++i) {
		wet_val += texture(wet_noise_tex, 0.04*freq*vpos).r/freq;
		freq    *= 2.0;
	}
	return sqrt(min(1.0, 8.0*wetness))*min(1.0, pow((wetness + max(wet_val, 0.6) - 0.6), 8.0));
}

float get_reflect_weight(in vec3 view_dir, in vec3 ws_normal, in float reflectivity2, in float refract_ix) {
	return reflectivity2 * mix(get_fresnel_reflection(view_dir, ws_normal, 1.0, refract_ix), 1.0, metalness);
}

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
	if (enable_clip_plane_z && vpos.z < clip_plane_z) discard;
#ifdef TRIPLANAR_TEXTURE
	vec4 texel = lookup_triplanar_texture(vpos, normal, tex0, tex0, tex0);
#else
#ifdef ENABLE_PARALLAX_MAP
	vec4 texel = texture(tex0, apply_parallax_map()); // FIXME: tex coord offset should apply to normal maps as well
#else
	vec4 texel = texture(tex0, tc);
	//texel = mix(vec4(1.0), texel, 0.001); // for debugging untextured
#endif // ENABLE_PARALLAX_MAP
#endif // TRIPLANAR_TEXTURE

#ifdef ENABLE_GAMMA_CORRECTION
	texel.rgb = pow(texel.rgb, vec3(2.2)); // gamma correction
#endif // ENABLE_GAMMA_CORRECTION

#ifdef TEXTURE_ALPHA_MASK
	if (texel.a < 0.99) discard;
#endif

#ifndef NO_ALPHA_TEST
	if (texel.a < min_alpha) discard; // Note: assumes walls don't have textures with alpha < 1
#endif
	float alpha = gl_Color.a;

	if (do_cube_lt_atten || do_sphere_lt_atten) { // && light_atten > 0.0
		vec3 v_inc = normalize(camera_pos - vpos);
		float dist = 0.0;

		if (do_cube_lt_atten) { // account for light attenuating/reflecting semi-transparent materials
			vec3 far_pt = vpos - 100.0*v_inc; // move it far away
			pt_pair res = clip_line(vpos, far_pt, cube_bb);
			dist        = length(res.v1 - res.v2);
		}
		if (do_sphere_lt_atten) { // alpha is calculated from distance between sphere intersection points
			float wpdist = distance(vpos, sphere_center);
			float dp     = dot(v_inc, (vpos - sphere_center));
			float adsq   = dp*dp - wpdist*wpdist + sphere_radius*sphere_radius;
			dist         = sqrt(max(adsq, 0.0)); // should always intersect
		}
		alpha += (1.0 - alpha)*(1.0 - exp(-light_atten*dist));

		if (refract_ix != 1.0 && dot(normal, v_inc) > 0.0) { // entering ray in front surface
			float reflect_w = get_fresnel_reflection(v_inc, normalize(normal), 1.0, refract_ix);
			alpha = reflect_w + alpha*(1.0 - reflect_w); // don't have a reflection color/texture, so just modify alpha
		} // else exiting ray in back surface - ignore for now since we don't refract the ray
	}
#ifndef NO_ALPHA_TEST
	if (keep_alpha && (texel.a * alpha) <= min_alpha) discard;
#endif

	float wetness = wet_effect;
	if (use_water_coverage && wetness > 0.0) {wetness *= get_water_snow_coverage();} // doesn't look quite right
	float reflectivity2 = reflectivity;

#ifdef ENABLE_PUDDLES
	if (wetness > 0.0 && wetness < 1.0 && normal.z > 0.5) { // create puddles for partially wet top surfaces
		wetness       = get_puddle_val(wetness);
		reflectivity2 = wetness;
	}
#endif // ENABLE_PUDDLES

#ifdef USE_WINDING_RULE_FOR_NORMAL
	float normal_sign = winding_normal_sign*((!two_sided_lighting || gl_FrontFacing) ? 1.0 : -1.0); // two-sided lighting
#else
	float normal_sign = ((!two_sided_lighting || (dot(eye_norm, epos.xyz) < 0.0)) ? 1.0 : -1.0); // two-sided lighting
#endif
	float wet_surf_val = wetness*max(normal.z, 0.0); // only +z surfaces are wet; doesn't apply to spec shininess though
	vec4 base_color = mix(1.0, 0.25, wet_surf_val)*gl_Color;

#ifdef ENABLE_SNOW_COVERAGE
	if (snow_cov_amt > 0.0 && normal.z > 0.4) {
		// add in snow on top of water/texture, using ratio of lit color from base color to pick up lighting
		float snow_amt = snow_cov_amt*get_water_snow_coverage()*min(1.0, 6.0*(normal.z-0.4));
		base_color     = mix(base_color, vec4(1.0), snow_amt);
		texel          = mix(texel, vec4(0.9, 0.9, 1.0, 1.0), snow_amt);
		alpha          = mix(alpha, 1.0, snow_amt);
	}
#endif
	vec3 lit_color  = emission.rgb + emissive_scale*gl_Color.rgb;
	lit_color      += base_color.rgb * get_indir_lighting(normal_sign) * mix(1.0, 0.7, wet_surf_val);
	vec3 normal_s   = normal_sign*normal;

#ifdef ENABLE_GAMMA_CORRECTION
	lit_color.rgb = pow(lit_color.rgb, vec3(2.2)); // gamma correction
#endif // ENABLE_GAMMA_CORRECTION

	if (direct_lighting) { // directional light sources with no attenuation
		vec3 n = normalize(normal_sign*eye_norm);
		if (enable_light0) {lit_color += add_light0(n, normal_sign, base_color);} // sun
		if (enable_light1) {lit_color += add_light1(n, normal_sign, base_color);} // moon
		if (enable_light2) {ADD_LIGHT(2);} // lightning
	}
	if (enable_dlights) {add_dlights_bm_scaled(lit_color, vpos, normalize(normal_s), base_color.rgb, 1.0, normal_sign, wet_surf_val);} // dynamic lighting
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));

#ifdef ENABLE_GAMMA_CORRECTION
	color.rgb = pow(color.rgb, vec3(0.45)); // gamma correction
#endif // ENABLE_GAMMA_CORRECTION

#ifdef ENABLE_REFLECTIONS // should this be before or after multiplication with texel?
	if (normal_s.z > 0.5 && vpos.z < reflect_plane_ztop && vpos.z > reflect_plane_zbot) { // top surface
		vec3 ws_normal   = normalize(normal_s);
		float ripple_mag = wetness * clamp(2.0*(1.0 - 0.5*length(epos.xyz)), 0.0, 1.0);

		if (ripple_mag > 0.0 && rain_intensity > 0.0) {
			// Note: since reflections are only enabled on vertical surfaces with normals in +z (world space), and the ripple normal map defaults to +z, no transforms are necessary
			ws_normal = normalize(mix(ws_normal, get_ripple_normal(1.0*tc, 0.2*ripple_time, wetness*rain_intensity), ripple_mag));
		}
#ifdef USE_BUMP_MAP
		ws_normal = normalize(mix(get_bump_map_normal(), ws_normal, 0.5*wetness));
#endif
		// Note: this doesn't work for refact_ix == 1, so we choose an arbitrary value of 1.3 (metals are lower, dielectrics are higher)
		float reflect_w = get_reflect_weight(normalize(camera_pos - vpos), ws_normal, reflectivity2, ((refract_ix == 1.0) ? 1.3 : refract_ix)); // default is water
		vec4 proj_pos   = fg_ProjectionMatrix * epos;
		vec2 ref_tex_st = clamp(0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5), 0.0, 1.0);
		color.rgb       = mix(color.rgb, texture(reflection_tex, ref_tex_st).rgb*get_wet_specular_color(wetness), reflect_w);
	}
#endif // ENABLE_REFLECTIONS

#ifdef ENABLE_CUBE_MAP_REFLECT
#ifdef USE_BUMP_MAP
	vec3 ws_normal = get_bump_map_normal();
#else
	vec3 ws_normal = normalize(normal_s);
#endif // USE_BUMP_MAP
	vec3 view_dir   = normalize(camera_pos - vpos);
	float reflect_w = get_reflect_weight(view_dir, ws_normal, reflectivity2, ((refract_ix == 1.0) ? 1.5 : refract_ix)); // default is not water
	vec3 rel_pos    = vpos - cube_map_center;
	rel_pos         = max(vec3(-cube_map_near_clip), min(vec3(cube_map_near_clip), rel_pos)); // clamp to cube bounds
	vec3 ref_dir    = rel_pos + cube_map_near_clip*ws_normal; // position offset within cube (approx.)
	color.rgb       = mix(color.rgb, texture(reflection_tex, ref_dir).rgb*specular_color.rgb, reflect_w);
#endif // ENABLE_CUBE_MAP_REFLECT

#ifdef APPLY_BURN_MASK
	color = apply_black_body_burn_mask(color, tc);
#endif

#ifndef SMOKE_ENABLED
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif // NO_ALPHA_TEST
#ifndef NO_FOG
	vec4 fcolor = fog_color;
	
	if (indir_lighting) {
		vec3 scene_urc  = scene_llc + scene_scale;
		float scene_bb[6]; scene_bb[0]=scene_llc.x; scene_bb[1]=scene_urc.x; scene_bb[2]=scene_llc.y; scene_bb[3]=scene_urc.y; scene_bb[4]=scene_llc.z; scene_bb[5]=scene_urc.z;
		float view_dist = distance(vpos, camera_pos);
		vec3 end_pos    = camera_pos + (vpos - camera_pos)*(min(fog_end, view_dist)/view_dist);
		pt_pair cres    = clip_line(end_pos, camera_pos, scene_bb);
		float scene_len = distance(cres.v2, cres.v1)/distance(end_pos, camera_pos);
		float pixel_lum = get_luminance(indir_lookup(cres.v1, normal)); // camera pos
		if (!underwater) {pixel_lum = mix(pixel_lum, get_luminance(indir_lookup(cres.v2, normal)), 0.75);} // FIXME: use multiple steps?
		//pixel_lum = get_luminance(lit_color.rgb)/max(0.01, get_luminance(gl_Color.rgb));
		fcolor.rgb *= mix(1.0, min(2.0*pixel_lum, 1.0), scene_len);
	}
	// FIXME: more physically correct to clip the view ray by the distance traveled through the water,
	// but not all shaders use this flow (leaves, plants, scenery, etc.)
	float fog_length = (underwater ? get_underwater_fog_length(water_depth, camera_pos, vpos) : length(epos.xyz));
	vec4 fog_out     = apply_fog_ffc(color, fog_length*get_custom_fog_scale_epos(epos), fcolor); // apply standard fog
	color = (keep_alpha ? vec4(fog_out.rgb, color.a) : fog_out);
#endif // NO_FOG
#else // SMOKE_ENABLED
	pt_pair res = clip_line(vpos, camera_pos, smoke_bb);
	if (res.v1 != res.v2) {add_smoke_contrib(res.v1, res.v2, color);}
#ifndef NO_ALPHA_TEST
	//color.a = min(color.a, texel.a); // Note: assumes walls don't have textures with alpha < 1
	if (color.a <= min_alpha) discard;
#endif
#endif // SMOKE_ENABLED
	//color = vec4(pow(color.rgb, vec3(0.45)), color.a); // gamma correction, doesn't really look right, and should un-gamma-correct textures when loading
	fg_FragColor = color;
	//fg_FragColor.r = wet_effect; fg_FragColor.b = 1.0-wet_effect;

#ifdef USE_DEPTH_TRANSPARENCY
	float d_delta = log_to_linear_depth(get_depth_at_fragment()) - log_to_linear_depth(gl_FragCoord.z);
	fg_FragColor.a *= clamp(d_delta/depth_trans_bias, 0.0, 1.0);
#endif
}
