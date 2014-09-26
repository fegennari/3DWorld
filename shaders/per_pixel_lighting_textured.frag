uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform vec4  emission  = vec4(0,0,0,1);

in vec3 dlpos, dl_normal; // world space
in vec3 normal;
//in vec2 tc; // comes from detail_normal_map.part.frag

vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
	return apply_bump_map(light_dir, eye_pos, normal, 1.0);
}

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;
	
	vec4 epos    = fg_ModelViewMatrix * vec4(dlpos, 1.0);
	vec3 normal2 = (no_normalize ? normal : normalize(normal)); // renormalize
	vec3 color   = emission.rgb;
	if (enable_dlights) color += add_dlights(dlpos, normalize(dl_normal), vec3(1.0)).rgb; // dynamic lighting
	if (enable_light0 ) color += add_light_comp_pos_smap_light0(normal2, epos).rgb;
	if (enable_light1 ) color += add_light_comp_pos_smap_light1(normal2, epos).rgb;
	fg_FragColor = vec4(texel.rgb * color, texel.a * gl_Color.a); // use gl_Color alpha directly
#ifndef NO_FOG
	fg_FragColor = apply_fog_epos(fg_FragColor, epos);
#endif
}
