uniform float half_dxy;
uniform float max_indir_light     = 1.0;
uniform float indir_vert_offset   = 0.25;
uniform float hemi_lighting_scale = 0.5;
uniform float hemi_lighting_normal_scale = 1.0;
// Note: these next two come from dynamic_lighting.part, which must always be included first
//uniform vec3 scene_llc, scene_scale; // scene bounds (world space)
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color     = vec3(0.0);
uniform vec3 out_range_indir_color = vec3(0.1);

// HEMI_LIGHTING uniforms
uniform vec3 sky_color;
uniform sampler2D ground_tex;

#ifdef USE_ALT_SCENE_BOUNDS
uniform vec3 alt_scene_llc, alt_scene_scale; // scene bounds (world space)
#endif

in vec3 vpos, normal; // world space

// unused, but could be used for indir lighting with 6 directional components
float cube_light_lookup(in vec3 normal_in, in vec3 plus_xyz, in vec3 minus_xyz) {
	vec3 n = normalize(normal_in);
	n      = sign(n) * n * n;
	return dot(max(n, 0.0), plus_xyz) + dot(max(-n, 0.0), minus_xyz);
}

vec3 map_to_indir_space(in vec3 pos) {
#ifdef USE_ALT_SCENE_BOUNDS
	vec3 ret = (pos - alt_scene_llc)/alt_scene_scale;
#else
	vec3 ret = (pos - scene_llc)/scene_scale;
#endif
#ifdef ENABLE_OUTSIDE_INDIR_RANGE
	return ret; // no clamp
#else
	return clamp(ret, 0.0, 1.0); // should be in [0.0, 1.0] range
#endif
}

#ifdef DIRECTIONAL_INDIR_LIGHTING
uniform sampler3D dir_light_pos_tex, dir_light_neg_tex;

vec3 indir_lookup(in vec3 pos, in vec3 n) {
	vec3 spos  = map_to_indir_space(pos);
	vec3 color = texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
	return color * cube_light_lookup(n, texture(dir_light_pos_tex, spos.zxy).xyz, texture(dir_light_neg_tex, spos.zxy).xyz);
}
#else
vec3 indir_lookup(in vec3 pos, in vec3 n) {
	vec3 spos = map_to_indir_space(pos);
#ifdef ENABLE_OUTSIDE_INDIR_RANGE
	if (spos.x < -0.01 || spos.y < -0.01 || spos.z < -0.01 || spos.x > 1.01 || spos.y > 1.01 || spos.z > 1.01) {return out_range_indir_color;} // outside the volume
#endif
	return texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
}
#endif

vec3 get_indir_lighting(in float normal_sign) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting || hemi_lighting) {
#ifdef USE_BUMP_MAP_INDIR // USE_BUMP_MAP must also be set
		// convert tangent space to eye space; should be okay to use transpose() of TBN rather than inverse()
		vec3 n_eye = transpose(get_tbn(1.0, eye_norm)) * get_bump_map_normal();
		vec3 n = normalize(inverse(fg_NormalMatrix) * n_eye); // convert eye space to world space
		n = mix(normal, n, bump_map_mag);
#else
		vec3 n = normal;
#endif
		vec3 spos = vpos + (normal_sign*indir_vert_offset*half_dxy)*n; // move slightly away from the vertex

		if (hemi_lighting && hemi_lighting_scale > 0.0) {
			float sky_lum   = max(get_luminance(sky_color), 0.33); // or use indir_color?
			spos            = map_to_indir_space(spos);
			vec3 gnd_color  = sky_lum*textureLod(ground_tex, spos.xy, 4).rgb; // use 4th LOD mipmap (64x64)
			vec3 hemi_color = mix(gnd_color, sky_color, (0.5*normal_sign*hemi_lighting_normal_scale*normal.z + 0.5));
			indir_color     = mix(indir_color, hemi_color, hemi_lighting_scale); // blend between the two
		}
		if (indir_lighting) {
			vec3 indir_light = indir_lookup(spos, n);
			//indir_light    = pow(indir_light, vec3(0.45)); // gamma correction - does this work with negative lights?
			indir_color     += indir_light; // indirect lighting
#if 0 // add faked indir specular
			vec3 eye_vect   = normalize(fg_ModelViewMatrixInverse[3].xyz - vpos);
			vec3 eye_ref    = reflect(eye_vect, normal_sign*normal);
			vec3 spos_spec  = vpos + (indir_vert_offset*half_dxy)*eye_ref;
			float spec_mag  = 1.0*pow(clamp((1.0 - dot(eye_vect, normal)), 0.0, 1.0), 4.0);
			indir_color    += spec_mag*specular_color.rgb*indir_lookup(spos_spec, n);
#endif
		}
	}
	return min(indir_color, max_indir_light);
}

void add_indir_lighting(inout vec3 lit_color, in float normal_sign) {
	lit_color += gl_Color.rgb * get_indir_lighting(normal_sign);
}
