uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta, half_dxy;
uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform vec3 smoke_color;
uniform float light_atten = 0.0, refract_ix = 1.0;
uniform float cube_bb[6];

// clipped eye position, clipped vertex position
varying vec3 eye, vpos, normal, lpos0, vposl; // world space
// epos, eye_norm, and tex_coord come from bump_map.frag

const float SMOKE_SCALE = 0.25;

// Note: dynamic point lights use reflection vector for specular, and specular doesn't move when the eye rotates
//       global directional lights use half vector for specular, which seems to be const per pixel, and specular doesn't move when the eye translates
#define ADD_LIGHT(i) lit_color += add_pt_light_comp(n, epos, i).rgb

vec3 add_light0(in vec3 n) {
	vec3 eye_pos   = epos.xyz;
	vec3 light_dir = normalize(gl_LightSource[0].position.xyz - eye_pos); // Note: could drop the -eye_pos for a directional light
	float nscale   = (use_shadow_map ? ((dot(n, light_dir) > 0.0) ? get_shadow_map_weight_light0(epos, n) : 0.0) : 1.0); // back-facing test

#ifdef DYNAMIC_SMOKE_SHADOWS
	if (lpos0 != vposl && nscale > 0.0) {
		vec3 dir      = vposl - lpos0;
		vec3 pos      = (lpos0 - scene_llc)/scene_scale;
		vec3 delta    = normalize(dir)*step_delta/scene_scale;
		float nsteps  = length(dir)/step_delta;
		int num_steps = 1 + min(100, int(nsteps)); // round up
		float step_weight  = fract(nsteps);
		float smoke_sscale = SMOKE_SCALE*step_delta/half_dxy;
	
		// smoke volume iteration using 3D texture, light0 to vposl
		for (int i = 0; i < num_steps; ++i) {
			float smoke = smoke_sscale*texture3D(smoke_and_indir_tex, pos.zxy).a*step_weight;
			nscale     *= (1.0 - smoke);
			pos        += delta*step_weight; // should be in [0.0, 1.0] range
			step_weight = 1.0;
		}
	}
#endif
	return add_light_comp_pos0(nscale*n, epos).rgb;
}

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
#ifdef ENABLE_PARALLAX_MAP
	vec4 texel  = texture2D(tex0, apply_parallax_map()); // FIXME: tex coord offset should apply to normal maps as well
#else
	vec4 texel  = texture2D(tex0, tex_coord);
#endif
	//texel.rgb = pow(texel.rgb, vec3(2.2)); // gamma correction

#ifdef TEXTURE_ALPHA_MASK
	if (texel.a < 0.99) discard;
#endif

#ifndef NO_ALPHA_TEST
	if (texel.a < min_alpha) discard; // Note: assumes walls don't have textures with alpha < 1
#endif
	float alpha = gl_Color.a;

	if (do_lt_atten) {
		vec3 v_inc = normalize(eye - vpos);

		//if (light_atten > 0.0) { // account for light attenuating/reflecting semi-transparent materials
			vec3 far_pt   = vpos - 100.0*v_inc; // move it far away
			pt_pair res   = clip_line(vpos, far_pt, cube_bb);
			float dist    = length(res.v1 - res.v2);
			float atten   = exp(-light_atten*dist);
			alpha += (1.0 - alpha)*(1.0 - atten);
		//}
		if (refract_ix != 1.0 && dot(normal, v_inc) > 0.0) { // entering ray in front surface
			float reflect_w = get_fresnel_reflection(v_inc, normalize(normal), 1.0, refract_ix);
			alpha = reflect_w + alpha*(1.0 - reflect_w); // don't have a reflection color/texture, so just modify alpha
		} // else exiting ray in back surface - ignore for now since we don't refract the ray
	}
#ifndef NO_ALPHA_TEST
	if (keep_alpha && (texel.a * alpha) <= min_alpha) discard;
#endif

#ifdef USE_LIGHT_COLORS
	vec3 lit_color = vec3(0);
#else
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
#endif
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
		vec3 n = normalize(normal_sign*normal);
#ifdef USE_LIGHT_COLORS
		vec3 dlight_color = gl_Color.rgb;
#else
		vec3 dlight_color = gl_FrontMaterial.diffuse.rgb;
#endif
		lit_color += add_dlights(vpos, n, eye, dlight_color); // dynamic lighting
	}
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));
	//color.rgb = pow(color.rgb, vec3(0.45)); // gamma correction

#ifdef APPLY_BURN_MASK
	color = apply_burn_mask(color, tex_coord);
#endif

#ifndef SMOKE_ENABLED
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif
#ifndef NO_FOG
	color = apply_fog_epos(color, epos); // apply standard fog
#endif
#else
	pt_pair res = clip_line(vpos, eye, smoke_bb);
	vec3 eye_c  = res.v1;
	vec3 vpos_c = res.v2;
	
	if (eye_c != vpos_c) { // smoke code
		vec3 dir      = eye_c - vpos_c;
		vec3 pos      = (vpos_c - scene_llc)/scene_scale;
		float nsteps  = length(dir)/step_delta;
		vec3 delta    = dir/(nsteps*scene_scale);
		int num_steps = 1 + min(1000, int(nsteps)); // round up
		float step_weight  = fract(nsteps);
		float smoke_sscale = SMOKE_SCALE*step_delta/half_dxy;
	
		// smoke volume iteration using 3D texture, pos to eye
		for (int i = 0; i < num_steps; ++i) {
			// Note: we could also lookup dynamic lighting here, which would be very slow but would also look really nice
			vec4 tex_val = texture3D(smoke_and_indir_tex, pos.zxy); // rgba = {color.rgb, smoke}
			float smoke  = smoke_sscale*tex_val.a*step_weight;
			float alpha  = (keep_alpha ? color.a : ((color.a == 0.0) ? smoke : 1.0));
			float mval   = ((!keep_alpha && color.a == 0.0) ? 1.0 : smoke);
			color        = mix(color, vec4((tex_val.rgb * smoke_color), alpha), mval);
			pos         += delta*step_weight; // should be in [0.0, 1.0] range
			step_weight  = 1.0;
		}
	}
#ifndef NO_ALPHA_TEST
	//color.a = min(color.a, texel.a); // Note: assumes walls don't have textures with alpha < 1
	if (color.a <= min_alpha) discard;
#endif
#endif
	//color = vec4(pow(color.r, 0.45), pow(color.g, 0.45), pow(color.b, 0.45), color.a); // gamma correction, doesn't really look right
	gl_FragColor = color;
}
