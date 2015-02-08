uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta, step_delta_shadow;
uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float emissive_scale = 0.0;
uniform float smoke_const_add = 0.0;
uniform vec3 smoke_color, sphere_center;
uniform vec3 sun_pos; // used for dynamic smoke shadows line clipping
uniform float light_atten = 0.0, refract_ix = 1.0;
uniform float cube_bb[6], sphere_radius;
uniform vec4 emission = vec4(0,0,0,1);

//in vec3 vpos, normal; // world space, come from indir_lighting.part.frag
// epos, eye_norm, and tc come from bump_map.frag
// camera_pos comes from dynamic_lighting.part

const float SMOKE_SCALE = 0.25;

// Note: dynamic point lights use reflection vector for specular, and specular doesn't move when the eye rotates
//       global directional lights use half vector for specular, which seems to be const per pixel, and specular doesn't move when the eye translates
#define ADD_LIGHT(i) lit_color += add_pt_light_comp(n, epos, i).rgb

vec3 add_light0(in vec3 n) {
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
	return add_light_comp_pos0(nscale*n, epos).rgb;
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
			tex_val.rgb += add_dlights(dl_pos, norm_dir, vec3(1.0)); // normal points from vertex to eye, color is applied later
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
		float smoke  = smoke_sscale*(tex_val.a + smoke_const_add)*step_weight*sd_ratio;
		float alpha  = (keep_alpha ? color.a : ((color.a == 0.0) ? smoke : 1.0));
		float mval   = ((!keep_alpha && color.a == 0.0) ? 1.0 : smoke);
		color        = mix(color, vec4((tex_val.rgb * smoke_color), alpha), mval);
		pos         += delta*step_weight; // should be in [0.0, 1.0] range
		step_weight  = 1.0;
	} // for i
}

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
#ifdef ENABLE_PARALLAX_MAP
	vec4 texel  = texture(tex0, apply_parallax_map()); // FIXME: tex coord offset should apply to normal maps as well
#else
	vec4 texel  = texture(tex0, tc);
#endif
	//texel.rgb = pow(texel.rgb, vec3(2.2)); // gamma correction
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

	vec3 lit_color = emission.rgb + emissive_scale*gl_Color.rgb;
	add_indir_lighting(lit_color);
	//lit_color.rgb = pow(lit_color.rgb, vec3(2.2)); // gamma correction

#ifdef USE_WINDING_RULE_FOR_NORMAL
	float normal_sign = ((!two_sided_lighting || gl_FrontFacing) ? 1.0 : -1.0); // two-sided lighting
#else
	float normal_sign = ((!two_sided_lighting || (dot(eye_norm, epos.xyz) < 0.0)) ? 1.0 : -1.0); // two-sided lighting
#endif
	
	if (direct_lighting) { // directional light sources with no attenuation
		vec3 n = normalize(normal_sign*eye_norm);
		if (enable_light0) lit_color += add_light0(n);
		if (enable_light1) lit_color += add_light_comp_pos_smap_light1(n, epos).rgb;
		if (enable_light2) ADD_LIGHT(2);
		if (enable_light3) ADD_LIGHT(3);
		if (enable_light4) ADD_LIGHT(4);
		if (enable_light5) ADD_LIGHT(5);
		if (enable_light6) ADD_LIGHT(6);
		if (enable_light7) ADD_LIGHT(7);
	}
	if (enable_dlights) {
		lit_color += add_dlights(vpos, normalize(normal_sign*normal), gl_Color.rgb); // dynamic lighting
	}
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));
	//color.rgb = pow(color.rgb, vec3(0.45)); // gamma correction

#ifdef APPLY_BURN_MASK
	color = apply_burn_mask(color, tc);
#endif

#ifndef SMOKE_ENABLED
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif // NO_ALPHA_TEST
#ifndef NO_FOG
	vec4 fog_out;
	
	if (indir_lighting) {
		vec3 scene_urc    = scene_llc + scene_scale;
		float scene_bb[6]; scene_bb[0]=scene_llc.x; scene_bb[1]=scene_urc.x; scene_bb[2]=scene_llc.y; scene_bb[3]=scene_urc.y; scene_bb[4]=scene_llc.z; scene_bb[5]=scene_urc.z;
		float view_dist   = distance(vpos, camera_pos);
		vec3 end_pos      = camera_pos + (vpos - camera_pos)*(min(fog_end, view_dist)/view_dist);
		pt_pair cres      = clip_line(end_pos, camera_pos, scene_bb);
		float scene_len   = distance(cres.v2, cres.v1)/distance(end_pos, camera_pos);
		float pixel_lum;
		if (indir_lighting) {pixel_lum = mix(get_luminance(indir_lookup(cres.v1)), get_luminance(indir_lookup(cres.v2)), 0.75);}
		else {pixel_lum = get_luminance(lit_color.rgb)/max(0.01, get_luminance(gl_Color.rgb));}
		vec4 fcolor  = fog_color;
		fcolor.rgb  *= mix(1.0, min(2.0*pixel_lum, 1.0), scene_len);
		fog_out = apply_fog_ffc(color, length(epos.xyz)*get_custom_fog_scale_epos(epos), fcolor); // apply standard fog
	}
	else {
		fog_out = apply_fog_epos(color, epos); // apply standard fog
	}
	color = (keep_alpha ? vec4(fog_out.rgb, color.a) : fog_out);
#endif // NO_FOG
#else // SMOKE_ENABLED
#ifdef NO_CLIP_SMOKE
	add_smoke_contrib(vpos, camera_pos, color);
#else
	pt_pair res = clip_line(vpos, camera_pos, smoke_bb);
	if (res.v1 != res.v2) {add_smoke_contrib(res.v1, res.v2, color);}
#endif // NO_CLIP_SMOKE
#ifndef NO_ALPHA_TEST
	//color.a = min(color.a, texel.a); // Note: assumes walls don't have textures with alpha < 1
	if (color.a <= min_alpha) discard;
#endif
#endif // SMOKE_ENABLED
	//color = vec4(pow(color.r, 0.45), pow(color.g, 0.45), pow(color.b, 0.45), color.a); // gamma correction, doesn't really look right
	fg_FragColor = color;
}
