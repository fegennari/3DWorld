uniform float half_dxy;
uniform float indir_vert_offset = 0.25;
// Note: these next two come from dynamic_lighting.part, which must always be included first
//uniform vec3 scene_llc, scene_scale; // scene bounds (world space)
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0.0);

// HEMI_LIGHTING uniforms
uniform vec3 sky_color;
uniform sampler2D ground_tex;

in vec3 vpos, normal; // world space

// unused, but could be used for indir lighting with 6 directional components
float cube_light_lookup(in vec3 normal_in, in vec3 plus_xyz, in vec3 minus_xyz) {
	vec3 n = normalize(normal_in);
	n      = sign(n) * n * n;
	return dot(max(n, 0.0), plus_xyz) + dot(max(-n, 0.0), minus_xyz);
}

#ifdef DIRECTIONAL_INDIR_LIGHTING
uniform sampler3D dir_light_pos_tex, dir_light_neg_tex;

vec3 indir_lookup(in vec3 pos, in vec3 n) {
	vec3 spos  = clamp((pos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
	vec3 color = texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
	return color * cube_light_lookup(n, texture(dir_light_pos_tex, spos.zxy).xyz, texture(dir_light_neg_tex, spos.zxy).xyz);
}
#else
vec3 indir_lookup(in vec3 pos, in vec3 n) {
	vec3 spos = clamp((pos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
	return texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
}
#endif

vec3 get_indir_lighting(in float normal_sign) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting || hemi_lighting) {
#ifdef USE_BUMP_MAP_INDIR // USE_BUMP_MAP must also be set
		vec3 n_eye = inverse(get_tbn(1.0)) * get_bump_map_normal(); // convert tangent space to eye space
		vec3 n = normalize(inverse(fg_NormalMatrix) * n_eye); // convert eye space to world space
#else
		vec3 n = normal;
#endif
		vec3 spos = vpos + (normal_sign*indir_vert_offset*half_dxy)*n; // move slightly away from the vertex

		if (hemi_lighting) {
			float sky_lum   = max(get_luminance(sky_color), 0.33); // or use indir_color?
			spos            = clamp((spos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
			vec3 gnd_color  = sky_lum*textureLod(ground_tex, spos.xy, 4).rgb; // use 4th LOD mipmap (64x64)
			vec3 hemi_color = mix(gnd_color, sky_color, (0.5*normal_sign*normal.z + 0.5));
			indir_color     = mix(indir_color, hemi_color, 0.5); // blend between the two
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
	return indir_color;
}

void add_indir_lighting(inout vec3 lit_color, in float normal_sign) {
	lit_color += gl_Color.rgb * get_indir_lighting(normal_sign);
}
