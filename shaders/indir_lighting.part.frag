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

vec3 indir_lookup(in vec3 pos) {
	vec3 spos = clamp((pos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
	return texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
}

void add_indir_lighting(inout vec3 lit_color) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting || hemi_lighting) {
#ifdef USE_BUMP_MAP_INDIR // USE_BUMP_MAP must also be set
		vec3 n_eye = inverse(get_tbn(1.0)) * get_bump_map_normal(); // convert tangent space to eye space
		vec3 n = normalize(inverse(fg_NormalMatrix) * n_eye); // convert eye space to world space
#else
		vec3 n = normal;
#endif
		vec3 spos = vpos + (indir_vert_offset*half_dxy)*n; // move slightly away from the vertex

		if (hemi_lighting) {
			float sky_lum   = max(get_luminance(sky_color), 0.33); // or use indir_color?
			spos            = clamp((spos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
			vec3 gnd_color  = sky_lum*textureLod(ground_tex, spos.xy, 4).rgb; // use 4th LOD mipmap (64x64)
			vec3 hemi_color = mix(gnd_color, sky_color, (0.5*normal.z + 0.5));
			indir_color     = mix(indir_color, hemi_color, 0.5); // blend between the two
		}
		if (indir_lighting) {
			vec3 indir_light = indir_lookup(spos);
			//indir_light    = pow(indir_light, vec3(0.45)); // gamma correction
			indir_color     += indir_light; // indirect lighting
		}
	}
	lit_color += gl_Color.rgb * indir_color;
}
