uniform float tscale = 1.0;
uniform float crater_scale = 1.0;
uniform sampler2D tex0;

in vec3 world_space_pos;
//in vec3 vpos, eye_norm, normal; // come from bump_map.part.frag and triplanar_bump_map.part.frag
//in vec4 epos; // comes from bump_map.part.frag

void main()
{
	vec4 texel = lookup_triplanar_texture(tscale*vpos, normalize(normal), tex0, tex0, tex0);
	vec3 norm_normal = normalize(eye_norm);

#ifdef HAS_CRATERS
	if (crater_scale > 0.0 && dot(norm_normal, normalize(fg_LightSource[0].position.xyz - epos.xyz)) > 0.0) { // facing the sun (smoother transition?)
		adjust_normal_for_craters(norm_normal, vpos); // add craters by modifying the normal
	}
#endif
	vec3 color = vec3(0.0);
	color += calc_shadow_atten(world_space_pos)*add_pt_light_comp(norm_normal, epos, 0).rgb; // sun_diffuse
	color += add_pt_light_comp(norm_normal, epos, 1).rgb; // galaxy_ambient
	fg_FragColor = vec4(texel.rgb * clamp(color, 0.0, 1.0), texel.a * gl_Color.a); // use gl_Color alpha directly
}

