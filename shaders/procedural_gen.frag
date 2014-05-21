uniform float min_alpha      = 0.0;
uniform float bump_tex_scale = 1.0;

//varying vec3 vpos, normal; // world space, come from indir_lighting.part.frag
//varying vec3 eye_norm; // comes from bump_map.frag

vec4 get_epos() {return fg_ModelViewMatrix * vec4(vpos, 1.0);}

#ifdef USE_BUMP_MAP
vec3 get_bump_map_normal() {
	return normalize(lookup_triplanar_texture(bump_tex_scale*vpos, normalize(normal), bump_map, bump_map, bump_map).xyz * 2.0 - 1.0);
}
vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
	vec3 binormal = normalize(cross(vec3(1,0,0), eye_norm));
	vec3 tangent  = normalize(cross(eye_norm, binormal));
	mat3 TBN  = transpose(mat3(tangent, binormal, eye_norm));
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return get_bump_map_normal();
}
#endif // USE_BUMP_MAP

void main()
{
	vec3 norm_normal = normalize(normal);
	vec4 texel;

	if (z_top_test) {
		if (use_noise_tex) {
			texel = get_texture_val_z_test(norm_normal, vpos);
		}
		else {
			lookup_triplanar_texture_2sz(vpos, norm_normal, tex0, tex0, tex0, tex1) * color0;
		}
	}
	else {
		if (use_noise_tex) {
			texel = get_texture_val(norm_normal, vpos);
		}
		else {
			texel = lookup_triplanar_texture(vpos, norm_normal, tex0, tex0, tex0) * color0;
		}
	}
	float alpha = gl_Color.a;
	vec3 lit_color = vec3(0.0);
	add_indir_lighting(lit_color);
	vec4 epos = get_epos();
	
	// directional light sources with no attenuation (Note: could add other lights later)
	if (enable_light0)  lit_color += add_light_comp_pos_smap_light0(normalize(eye_norm), epos).rgb;
	if (enable_light1)  lit_color += add_light_comp_pos_smap_light1(normalize(eye_norm), epos).rgb;
	if (enable_dlights) lit_color += add_dlights(vpos, norm_normal, gl_Color.rgb); // dynamic lighting
	fg_FragColor = vec4((texel.rgb * lit_color), (texel.a * alpha));
#ifndef NO_ALPHA_TEST
	if (fg_FragColor.a <= min_alpha) discard;
#endif
#ifndef NO_FOG
	fg_FragColor = apply_fog_epos(fg_FragColor, epos); // apply standard fog
#endif
}
