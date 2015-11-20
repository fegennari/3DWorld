uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform float opacity   = 1.0;
uniform vec4 color_scale = vec4(1.0);
uniform float smap_atten_cutoff = 10.0; // for TT mode
uniform float smap_atten_slope  = 0.5;  // for TT mode

in vec2 tc;
in vec4 epos;
in vec3 normal; // eye space
in vec3 ws_pos;
in vec3 ws_normal;

// for indir lighting
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0.0);

void add_indir_lighting(inout vec3 lit_color, in vec3 vpos) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting) {
		vec3 spos    = clamp((vpos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		indir_color += texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
	}
	lit_color += gl_Color.rgb * indir_color;
}

void main() {
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;
#ifdef ENABLE_OPACITY
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
#endif
	vec3 color = vec3(0.0);
	float ambient_scale = (indir_lighting ? 0.0 : 1.0);
	float smap_scale = 1.0;
#ifdef USE_SMAP_SCALE
	smap_scale = clamp(smap_atten_slope*(smap_atten_cutoff - length(epos.xyz)), 0.0, 1.0);
#endif
	add_indir_lighting(color, ws_pos);
	if (enable_light0)  {color += add_leaf_light_comp(normal, epos, 0, ambient_scale, ((use_shadow_map && smap_scale > 0.0) ? mix(1.0, get_shadow_map_weight_light0(epos, normal), smap_scale) : 1.0)).rgb;}
	if (enable_light1)  {color += add_leaf_light_comp(normal, epos, 1, ambient_scale, ((use_shadow_map && smap_scale > 0.0) ? mix(1.0, get_shadow_map_weight_light1(epos, normal), smap_scale) : 1.0)).rgb;}
	if (enable_light2)  {color += add_pt_light_comp  (normalize(normal), epos, 2).rgb;} // lightning
	if (enable_dlights) {add_dlights(color, ws_pos, ((dot(normal, epos.xyz) > 0.0) ? -1.0 : 1.0)*normalize(ws_normal), vec3(1.0));}
	fg_FragColor = texel*vec4(min(2.0*gl_Color.rgb, clamp(color*color_scale.rgb, 0.0, 1.0)), 1.0); // limit lightning color
#ifndef NO_FOG
	fg_FragColor = apply_fog(fg_FragColor);
#endif
}
