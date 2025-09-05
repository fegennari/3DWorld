#include <fresnel.part>

#ifdef USE_VOROCRACKS
#include <vorocracks.part.frag>
#endif

uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D wet_noise_tex;
uniform sampler2D sky_zval_tex;
uniform sampler2D emissive_map; // only used when ENABLE_EMISSIVE_MAP
uniform sampler2D blue_noise_tex; // used for smoke shadow maps
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
uniform float diffuse_scale = 1.0; // also applies to specular
uniform float cube_bb[6], sphere_radius;
uniform float depth_trans_bias, clip_plane_z, ripple_time, rain_intensity, reflectivity, snow_cov_amt;
uniform float reflect_plane_ztop, reflect_plane_zbot;
uniform float winding_normal_sign = 1.0;
uniform float cube_map_reflect_mipmap_level = 0.0;
uniform float puddle_scale        = 0.04;
uniform float water_damage_zmax   = 0.0;
uniform float water_damage_zscale = 1.0;
uniform float crack_scale  = 1.0;
uniform float crack_sharp  = 100.0;
uniform float crack_weight = 0.0;
uniform float crack_zmax   = 0.0;
uniform float crack_normal_zmax   = -2.0; // default is all normals
uniform int cube_map_texture_size = 0;
uniform vec4 emission = vec4(0,0,0,1);
uniform bool two_sided_lighting = true;

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
#ifdef ENABLE_SKY_OCCLUSION
uniform float sky_occlude_scale = 0.5;
#endif

const float SMOKE_SCALE = 0.25;

vec3 get_closest_cube_normal(in float cube[6], in vec3 pos) { // assumes pos has been clipped to cube; cube = {-x, +x, -y, +y, -z, +z}
#if 0 // smoother edges
	vec3 center = 0.5*vec3(cube[0]+cube[1], cube[2]+cube[3], cube[4]+cube[5]);
	return normalize(sign(pos - center)*pow(abs(normalize(pos - center)), vec3(80.0)));
#else
	float xn = abs(pos.x - cube[0]), xp = abs(pos.x - cube[1]);
	float yn = abs(pos.y - cube[2]), yp = abs(pos.y - cube[3]);
	float zn = abs(pos.z - cube[4]), zp = abs(pos.z - cube[5]);
	vec3 n   = vec3(-1,0,0);
	float d  = xn; // -x
	if (xp < d) {d = xp; n = vec3( 1,0,0);} // +x
	if (yn < d) {d = yn; n = vec3(0,-1,0);} // -y
	if (yp < d) {d = yp; n = vec3(0, 1,0);} // +y
	if (zn < d) {d = zn; n = vec3(0,0,-1);} // -z
	if (zp < d) {d = zp; n = vec3(0,0, 1);} // +z
	return n;
#endif
}

void ray_trace_cube_sphere(in vec3 p1, in vec3 dir, out vec3 p2, out vec3 n) {
	if (light_atten > 0.0) { // cube case
		vec3 far_pt = p1 + 100.0*dir; // move it far away
		pt_pair res = clip_line(p1, far_pt, cube_bb);
		p2 = res.v1;
		n  = get_closest_cube_normal(cube_bb, p2); // use tmax point
	}
	else { // sphere case
		vec3 ray_dir = normalize(p1 - sphere_center)*sphere_radius; // renormalize to sphere radius to reduce distortion
		float dist   = abs(dot(dir, ray_dir));
		p2 = (ray_dir + sphere_center) + dist*dir; // other intersection point
		n  = normalize(p2 - sphere_center);
	}
}

// Note: dynamic point lights use reflection vector for specular, and specular doesn't move when the eye rotates
//       global directional lights use half vector for specular, which seems to be const per pixel, and specular doesn't move when the eye translates
#define ADD_LIGHT(i) lit_color += add_pt_light_comp(n, epos, i).rgb

float get_ambient_scale() {
#ifdef ENABLE_SKY_OCCLUSION
	vec3 pos    = (vpos - scene_llc)/scene_scale;
	float cmp   = vpos.z + 0.5*half_dxy;
	float delta = texture(sky_zval_tex, pos.xy).r - cmp; // incorrectly interpolated, but smooth
	return mix(1.0, clamp((1.0 - 0.1*delta/half_dxy), 0.0, 1.0), sky_occlude_scale);
#else
	return 1.0;
#endif
}

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
		float nsteps  = length(dir)/step_delta;
		int num_steps = 1 + min(100, int(nsteps)); // round up
		float sd_ratio= max(1.0, nsteps/100.0);
		vec3 delta    = normalize(dir)*sd_ratio*step_delta/scene_scale;
		float step_weight  = fract(nsteps);
		float smoke_sscale = SMOKE_SCALE*sd_ratio*step_delta/half_dxy;
	
		// smoke volume iteration using 3D texture, light0 to vposl
		for (int i = 0; i < num_steps; ++i) {
			float smoke = smoke_sscale*texture(smoke_and_indir_tex, pos.zxy).a*step_weight;
			nscale     *= (1.0 - smoke);
			pos        += delta*step_weight; // should be in [0.0, 1.0] range
			step_weight = 1.0;
		}
	}
#endif // DYNAMIC_SMOKE_SHADOWS
	return add_light_comp_pos_scaled_light(nscale*n, epos, diffuse_scale, get_ambient_scale(), base_color, fg_LightSource[0], normal_sign).rgb;
}

vec3 add_light1(in vec3 n, in float normal_sign, in vec4 base_color) {
#ifdef USE_SHADOW_MAP
	if (use_shadow_map) {n *= get_shadow_map_weight_light1(epos, n);}
#endif
	return add_light_comp_pos_scaled_light(n, epos, diffuse_scale, get_ambient_scale(), base_color, fg_LightSource[1], normal_sign).rgb;
}

void add_smoke_contrib(in vec3 eye_c, in vec3 vpos_c, inout vec4 color) {
	vec3 dir      = eye_c - vpos_c;
	vec3 norm_dir = normalize(dir); // used for dlights
	vec3 pos      = (vpos_c - scene_llc)/scene_scale;
	float nsteps  = length(dir)/step_delta;
	int num_steps = 1 + min(256, int(nsteps)); // round up
	float sd_ratio= max(1.0, nsteps/256.0);
	vec3 delta    = sd_ratio*dir/(nsteps*scene_scale);
	float step_weight  = fract(nsteps);
	float smoke_sscale = SMOKE_SCALE*step_delta/half_dxy;
#if defined(SMOKE_SHADOW_MAP) && defined(USE_SHADOW_MAP)
	vec4 cur_epos   = fg_ModelViewMatrix * vec4(vpos_c, 1.0);
	vec3 epos_delta = fg_NormalMatrix * (delta * scene_scale); // eye space pos/delta
	cur_epos.xyz   += texture(blue_noise_tex, gl_FragCoord.xy/textureSize(blue_noise_tex, 0)).r*epos_delta; // add blue noise offset to break up banding artifacts
	const float smoke_albedo = 0.9;
	vec3 sl_color   = smoke_albedo * fg_LightSource[0].diffuse.rgb;
#endif

	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec4 tex_val = texture(smoke_and_indir_tex, pos.zxy); // rgba = {color.rgb, smoke}
#ifdef SMOKE_DLIGHTS
		if (enable_dlights) { // dynamic lighting
			vec3 dl_pos = pos*scene_scale + scene_llc;
			add_dlights_bm_scaled(tex_val.rgb, dl_pos, epos, norm_dir, vec3(1.0), 0.0, 1.0, 0.0); // normal points from vertex to eye, override bump mapping, color is applied later
		}
#endif // SMOKE_DLIGHTS
#if defined(SMOKE_SHADOW_MAP) && defined(USE_SHADOW_MAP)
		if (enable_light0) {
			tex_val.rgb  += sl_color * get_shadow_map_weight_light0_no_bias(cur_epos);
			cur_epos.xyz += epos_delta*step_weight;
		}
#endif // USE_SHADOW_MAP && SMOKE_SHADOW_MAP

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

float noise_lookup_4_octaves(in vec3 pos) {
	float val  = 0.0;
	float freq = 1.0;

	for (int i = 0; i < 4; ++i) {
		val  += texture(wet_noise_tex, freq*pos).r/freq;
		freq *= 2.0;
	}
	return val;
}

float get_puddle_val(in float wetness) {
	float wet_val = noise_lookup_4_octaves(puddle_scale*vec3(vpos.x, vpos.y, vpos.z*water_damage_zscale));
	return sqrt(min(1.0, 8.0*wetness))*min(1.0, pow((wetness + max(wet_val, 0.6) - 0.6), 8.0));
}

float get_crack_factor(in float weight) {
#ifdef USE_VOROCRACKS
	float val = get_crack_weight(crack_scale*tc);
	float sparse_noise_scale = 0.435; // lower frequency
#else
	float val = noise_lookup_4_octaves(crack_scale*(vpos + vec3(1.234))); // offset so that it's not aligned to puddles
	val = min(1.0, crack_sharp*crack_scale*abs(val - 1.0));
	float sparse_noise_scale = 1.237;
#endif
	float weight0 = 0.5;

	if (weight0 < 1.0) { // handle weight less than 1.0 with sparse cracks
		float strength = noise_lookup_4_octaves(sparse_noise_scale*crack_scale*(vpos + vec3(4.321)));
		float sub_val  = 1.875*(0.5*weight0 + 0.25); // 1.875 = sum of 4 weight octaves
		val = mix(val, 1.0, clamp(20.0*(strength - sub_val), 0.0, 1.0)); // sharp cutoff
	}
	return val;
}

// returns {fresnel_term, reflect_weight}
vec2 get_reflect_weight(in vec3 view_dir, in vec3 ws_normal, in float reflectivity2, in float refract_ix) {
	float fresnel = get_fresnel_reflection(view_dir, ws_normal, 1.0, refract_ix);
	return vec2(fresnel, (reflectivity2 * mix(fresnel, 1.0, metalness)));
}

#ifdef ENABLE_CUBE_MAP_REFLECT
vec3 cube_map_reflect_texture_filter(in vec3 ref_dir, in int mip_bias, in int filt_sz) {

	float level = max(mip_bias, textureQueryLod(reflection_tex, ref_dir).y);
	vec3 rdir   = normalize(ref_dir);
	vec3 d0     = vec3(0);

	if (abs(rdir.x) < abs(rdir.y)) {
		if (abs(rdir.x) < abs(rdir.z)) {d0.x = 1;} else {d0.z = 1;}
	}
	else {
		if (abs(rdir.y) < abs(rdir.z)) {d0.y = 1;} else {d0.z = 1;}
	}
	float exp_scale = 1.0/filt_sz;
	float texel_sz  = pow(2.0, mip_bias)/cube_map_texture_size; // size of one texel for current mip level
	vec3 d1         = texel_sz*normalize(cross(d0, rdir));
	vec3 d2         = texel_sz*normalize(cross(d1, rdir));
	float tot_w     = 0.0;
	vec3 color      = vec3(0.0);

	for (int y = -filt_sz; y <= filt_sz; ++y) {
		vec3 dir = (rdir - filt_sz*d1 + y*d2);
		for (int x = -filt_sz; x <= filt_sz; ++x) {
			float weight = exp(-exp_scale*sqrt(x*x + y*y)); // Gaussian
			color       += weight*textureLod(reflection_tex, dir, level).rgb; // no anisotropic filtering
			tot_w       += weight;
			dir         += d1;
		}
	}
	return color/tot_w;
}

vec3 apply_cube_map_blur(in vec3 ref_dir, in int blur_val) { // blur_val ranges from 0 to around 8
	if (blur_val <= 0) {return texture(reflection_tex, ref_dir).rgb;} // auto mipmaps + anisotropic filtering with no blur
	// add some blur
	int mip_bias = 0;
	int filt_sz  = 1;

	// split blur logarithmically across mip bias and filter kernel size, to a max filter size of 7 (15x15)
	// this provides the best quality vs. performance trade-off (mip bias is free but blocky and has seams; filter is smooth but expensive)
	while (blur_val > 0) {
		// Note: fliter_sz from 1=>2 is really a 1x1 to 3x3 transition, by we still use 2x to account for the Gaussian weights vs. box filter used for mipmaps
		if (filt_sz  < 8) {--blur_val; filt_sz *= 2;}
		if (blur_val > 0) {--blur_val; ++mip_bias;}
	}
	--filt_sz; // filt_sz of 0 is a 1x1 box
	if      (filt_sz == 1) {return cube_map_reflect_texture_filter(ref_dir, mip_bias, 1);}
	else if (filt_sz == 3) {return cube_map_reflect_texture_filter(ref_dir, mip_bias, 3);}
	else if (filt_sz == 5) {return cube_map_reflect_texture_filter(ref_dir, mip_bias, 5);}
	else                   {return cube_map_reflect_texture_filter(ref_dir, mip_bias, 7);}
}
#endif // ENABLE_CUBE_MAP_REFLECT

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main() {

	if (enable_clip_plane_z && vpos.z < clip_plane_z) discard;
#ifdef ENABLE_ABSTRACT_ART
	vec4 texel = gen_abstract_art(tc, gl_Color.rgb); // color RGB is the seed
#else
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
#endif // ENABLE_ABSTRACT_ART

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

	if (do_lt_atten && light_atten != 0.0) {
		vec3 v_inc = normalize(camera_pos - vpos);
		float dist = 0.0;

		if (light_atten > 0.0) { // cube case; account for light attenuating/reflecting semi-transparent materials
			vec3 far_pt = vpos - 100.0*v_inc; // move it far away
			pt_pair res = clip_line(vpos, far_pt, cube_bb);
			dist        = light_atten*length(res.v1 - res.v2);
		}
		else { // sphere case; alpha is calculated from distance between sphere intersection points
			dist = abs(dot(v_inc, normalize(vpos - sphere_center)*sphere_radius));
			dist = -light_atten*dist; // should always intersect; note that light_atten is negative
		}
		alpha += (1.0 - alpha)*(1.0 - exp(-dist));

#ifndef ENABLE_CUBE_MAP_REFLECT // skip refraction for cube maps since we don't want to reflect the internal surface medium
		if (refract_ix != 1.0 && dot(normal, v_inc) > 0.0) { // entering ray in front surface
			float reflect_w = get_fresnel_reflection(v_inc, normalize(normal), 1.0, refract_ix); // air => object transition
			alpha = reflect_w + alpha*(1.0 - reflect_w); // don't have a reflection color/texture, so just modify alpha
		} // else exiting ray in back surface - ignore for now since we don't refract the ray
#endif // not ENABLE_CUBE_MAP_REFLECT
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

#ifdef ENABLE_WATER_DAMAGE
	if (wetness > 0.0) { // add water damage similar to puddles; doesn't apply to emissive surfaces
		wetness       = ((vpos.z < water_damage_zmax && emissive_scale == 0.0) ? get_puddle_val(wetness) : 0.0);
		reflectivity2 = wetness;
	}
#endif // ENABLE_WATER_DAMAGE

#ifdef ADD_CRACKS
	if (crack_weight > 0.0 && vpos.z < crack_zmax && normal.z > crack_normal_zmax) {
		float color_scale = get_crack_factor(crack_weight);
		texel.rgb *= color_scale;
		wetness   *= color_scale;
	}
#endif // ADD_CRACKS

#ifdef USE_WINDING_RULE_FOR_NORMAL
	float normal_sign  = winding_normal_sign*((!two_sided_lighting || gl_FrontFacing) ? 1.0 : -1.0); // two-sided lighting
#else
	float normal_sign  = ((!two_sided_lighting || (dot(eye_norm, epos.xyz) < 0.0)) ? 1.0 : -1.0); // two-sided lighting
#endif
	vec3 normal_s      = normal_sign*normal;
#ifdef ENABLE_WATER_DAMAGE
	float wet_surf_val = wetness; // all surfaces are wet
#else
	float wet_surf_val = wetness*max(normal.z, 0.0); // only +z surfaces are wet; doesn't apply to spec shininess though
#endif
	
#ifdef ENABLE_ABSTRACT_ART
	vec4 base_color    = vec4(1.0); // white
#else
	vec4 base_color    = gl_Color; // unlit material color (albedo)
#endif // ENABLE_ABSTRACT_ART

	base_color        *= mix(1.0, 0.25, wet_surf_val); // apply wet surface darkening
	base_color.rgb    *= texel.rgb; // maybe gamma correction is wrong when the texture is added this way? but it's not actually used in any scenes
	vec3 lit_color     = emission.rgb*texel.rgb + emissive_scale*base_color.rgb;

#ifdef ENABLE_SNOW_COVERAGE
	if (snow_cov_amt > 0.0 && normal_s.z > 0.4) {
		// add in snow on top of water/texture, using ratio of lit color from base color to pick up lighting
		float snow_amt = snow_cov_amt*get_water_snow_coverage()*min(1.0, 6.0*(normal_s.z-0.4));
		base_color     = mix(base_color, vec4(1.0), snow_amt);
		texel          = mix(texel, vec4(0.9, 0.9, 1.0, 1.0), snow_amt);
		alpha          = mix(alpha, 1.0, snow_amt);
	}
#endif // ENABLE_SNOW_COVERAGE
#ifdef ENABLE_EMISSIVE_MAP
	lit_color     += texture(emissive_map, tc);
#endif
	lit_color     += base_color.rgb * get_indir_lighting(normal_sign) * mix(1.0, 0.7, wet_surf_val);

#ifdef ENABLE_GAMMA_CORRECTION
	lit_color.rgb = pow(lit_color.rgb, vec3(2.2)); // gamma correction
#endif // ENABLE_GAMMA_CORRECTION

	if (direct_lighting) { // directional light sources with no attenuation
		vec3 n = normalize(normal_sign*eye_norm);
		if (enable_light0) {lit_color += add_light0(n, normal_sign, base_color);} // sun
		if (enable_light1) {lit_color += add_light1(n, normal_sign, base_color);} // moon
		if (enable_light2) {ADD_LIGHT(2);} // lightning
	}
	if (enable_dlights) {add_dlights_bm_scaled(lit_color, vpos, epos, normalize(normal_s), base_color.rgb, 1.0, normal_sign, wet_surf_val);} // dynamic lighting
	vec4 color = vec4(lit_color, (texel.a * alpha));

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
		vec3 spec_scale = vec3(1.0);
#ifdef USE_SPEC_MAP
		spec_scale     *= get_spec_color();
#endif
		// Note: this doesn't work for refact_ix == 1, so we choose an arbitrary value of 1.3 (metals are lower, dielectrics are higher)
		float reflect_w = get_reflect_weight(normalize(camera_pos - vpos), ws_normal, reflectivity2, ((refract_ix == 1.0) ? 1.3 : refract_ix)).y; // default is water
		vec4 proj_pos   = fg_ProjectionMatrix * epos;
		vec2 ref_tex_st = clamp(0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5), 0.0, 1.0);
		color.rgb       = mix(color.rgb, texture(reflection_tex, ref_tex_st).rgb*get_wet_specular_color(wetness), reflect_w*spec_scale);
	}
#endif // ENABLE_REFLECTIONS

#ifdef ENABLE_CUBE_MAP_REFLECT
#ifdef ENABLE_BUILDING_CUBE_MAP // handle metal and glass reflections
	if ((metalness > 0.0 && specular_color.rgb != vec3(0.0) && gl_Color.rgb != vec3(0.0)) || refract_ix > 1.0) {
		vec3 view_dir  = normalize(vpos - camera_pos);
		vec3 ws_normal = normalize(normal_s);
		vec3 ref_dir   = reflect(view_dir, ws_normal);
		vec3 boxMin    = cube_map_center - cube_map_near_clip;
		vec3 boxMax    = cube_map_center + cube_map_near_clip;
		vec3 firstPlaneIntersect  = (boxMax - vpos) / ref_dir;
		vec3 secondPlaneIntersect = (boxMin - vpos) / ref_dir;
		vec3 furthestPlane = max(firstPlaneIntersect, secondPlaneIntersect);
		float dist         = min(min(furthestPlane.x, furthestPlane.y), furthestPlane.z);
		vec3 intersectPosition = vpos + ref_dir * dist;
		ref_dir = intersectPosition - cube_map_center;

		if (metalness > 0.0) { // metal
			float shininess= specular_color.a; // typically 1-80
			int blur_val   = max(0, int((80.0 - shininess)/10.0));
			vec3 ref_tex   = apply_cube_map_blur(ref_dir, blur_val);
			// white specular: modulate with material color (for different shades of metal)
			vec3 spec_color= ((specular_color.r == specular_color.g && specular_color.r == specular_color.b) ? texel.rgb*gl_Color.rgb : specular_color.rgb);
			color.rgb = mix(color.rgb, spec_color*ref_tex, metalness);
		}
		else { // glass/dielectric
			vec2 reflected = get_reflect_weight(-view_dir, ws_normal, reflectivity2, refract_ix); // {fresnel_term, reflect_weight}
			color = mix(color, vec4(texture(reflection_tex, ref_dir).rgb, 1.0), reflected.y);
		}
	}
#else // !ENABLE_BUILDING_CUBE_MAP
#ifdef ENABLE_CUBE_MAP_BUMP_MAPS
	vec3 ws_normal = get_bump_map_normal();
#else
	vec3 ws_normal = normalize(normal_s);
#endif // ENABLE_CUBE_MAP_BUMP_MAPS
	vec3 spec_scale = vec3(1.0);
#ifdef USE_SPEC_MAP
	spec_scale     *= get_spec_color();
#endif
	float ref_ix    = refract_ix;
	vec3 view_dir   = normalize(vpos - camera_pos);
	vec2 reflected  = get_reflect_weight(-view_dir, ws_normal, reflectivity2, ref_ix); // {fresnel_term, reflect_weight}
	vec3 reflect_w  = reflected.y*spec_scale; // use reflected weight including metalness
	vec3 rel_pos    = vpos - cube_map_center;
	rel_pos         = max(vec3(-cube_map_near_clip), min(vec3(cube_map_near_clip), rel_pos)); // clamp to cube bounds
	vec3 ref_dir    = rel_pos + cube_map_near_clip*reflect(view_dir, ws_normal); // position offset within cube (approx.)
	//vec3 ref_dir    = view_dir; // invisible effect
	vec3 t_color    = color.rgb; // transmitted color
	
	if (alpha < 1.0) { // handle refraction
		vec3 refract_dir   = refract(view_dir, ws_normal, 1.0/ref_ix); // refraction into the object
		vec3 int_ref_color = vec3(0);
		float int_ref_w    = 0.0;

		if (do_lt_atten && light_atten != 0.0 && refract_ix != 1.0) {
			vec3 p2, ref_n;
			ray_trace_cube_sphere(vpos, refract_dir, p2, ref_n);

			if (light_atten < 0.0 || dot(ref_n, refract_dir) > 0.0) { // camera facing test for cube faces
				// update cube map lookup position to the back side where the refracted ray exits
				vec3 rel_pos2 = p2 - cube_map_center; // relative pos on back side
				vec3 ref_dir2 = reflect((rel_pos2 + cube_map_near_clip*refract_dir), -ref_n); // internal reflection vector
				// compute refraction parameters
				refract_dir = refract(refract_dir, -ref_n, ref_ix/1.0); // refraction out of the object
				int_ref_w   = get_fresnel_reflection(-refract_dir, -ref_n, refract_ix, 1.0); // object => air transition
				refract_dir = rel_pos2 + cube_map_near_clip*refract_dir; // refract dir normalized to cube map
				// compute internal reflection
				int_ref_color = texture(reflection_tex, ref_dir2).rgb;
#if 0
				vec3 p2b, ref_nb;
				ray_trace_cube_sphere(p2, ref_dir2, p2b, ref_nb);
				vec3 refract_dir2 = refract(ref_dir2, ref_nb, ref_ix/1.0); // refraction out of the object
#endif
			}
		}
		t_color = mix(mix(texture(reflection_tex, refract_dir).rgb, int_ref_color, int_ref_w), t_color, alpha); // not multiplied by specular color
		alpha   = 1.0;
	}
	int blur_val = int(cube_map_reflect_mipmap_level + 0.5); // ranges from 1.0 (no blur) to around 256 (max blur)
	color.rgb    = apply_cube_map_blur(ref_dir, blur_val);
	// use fresnel term to select between metallic specular color and monochrome fresnel specular color
	vec3 spec_color = mix(specular_color.rgb, vec3((specular_color.r + specular_color.g + specular_color.b)/3.0), reflected.x);
	color.rgb       = mix(t_color, color.rgb*spec_color, reflect_w);
#endif // ENABLE_BUILDING_CUBE_MAP
#endif // ENABLE_CUBE_MAP_REFLECT

#ifdef APPLY_BURN_MASK
	if (burn_offset > -1.0) {color = apply_black_body_burn_mask(color, tc);}
#endif

#ifndef SMOKE_ENABLED
#ifndef NO_ALPHA_TEST
	if (color.a <= min_alpha) discard;
#endif // NO_ALPHA_TEST
#ifndef NO_FOG
	vec4 fcolor = fog_color;
	
	if (indir_lighting) {
#ifdef USE_SIMPLE_INDIR
		vec3 pixel_color = indir_lookup(vpos, normal);
#else
		vec3 scene_urc  = scene_llc + scene_scale;
		float scene_bb[6]; scene_bb[0]=scene_llc.x; scene_bb[1]=scene_urc.x; scene_bb[2]=scene_llc.y; scene_bb[3]=scene_urc.y; scene_bb[4]=scene_llc.z; scene_bb[5]=scene_urc.z;
		float view_dist = distance(vpos, camera_pos);
		vec3 end_pos    = camera_pos + (vpos - camera_pos)*(min(fog_end, view_dist)/view_dist);
		pt_pair cres    = clip_line(end_pos, camera_pos, scene_bb);
		float scene_len = distance(cres.v2, cres.v1)/distance(end_pos, camera_pos);
		vec3 pixel_color;
		if (underwater) {pixel_color = indir_lookup(cres.v1, normal);} // camera pos
		else {
			pixel_color = indir_lookup(cres.v2, normal); // target pos
			const int num_steps = 4; // + starting step
			vec3 delta = (cres.v1 - cres.v2)/num_steps;
			vec3 pos = cres.v2;
			for (int i = 0; i < num_steps; ++i) {
				pos        += delta;
				pixel_color = mix(pixel_color, indir_lookup(pos, normal), 0.25);
			}
		}
		//pixel_color = lit_color.rgb/max(vec3(0.01), gl_Color.rgb);
		fcolor.rgb *= mix(vec3(1.0), min(2.0*pixel_color, vec3(1.0)), scene_len);
#endif // USE_SIMPLE_INDIR
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
	float d_delta   = log_to_linear_depth(get_depth_at_fragment()) - log_to_linear_depth(gl_FragCoord.z);
	fg_FragColor.a *= smoothstep(0.0, 1.0, clamp(d_delta/depth_trans_bias, 0.0, 1.0));
#endif
}
