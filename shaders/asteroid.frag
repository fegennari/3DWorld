varying vec3 vpos, normal, world_normal;
varying vec4 epos;

void main()
{
	vec3 n     = normalize(normal);
	vec4 texel = get_texture_val(normalize(world_normal), vpos);
	vec4 color = gl_FrontMaterial.emission;

	for (int i = 0; i < 8; ++i) {
		color += add_pt_light_comp(n, epos, i);
	}
	color = vec4(texel.rgb * clamp(color.rgb, 0.0, 1.0), texel.a * gl_Color.a); // use diffuse alpha directly
#ifndef NO_FOG
	color = apply_fog(color); // apply standard fog
#endif
	gl_FragColor = color;
}

