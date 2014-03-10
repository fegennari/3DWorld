uniform sampler3D ao_tex, shadow_tex;

varying vec3 vpos, normal, world_normal;
varying vec4 epos;

void main()
{
	vec3 n     = normalize(normal);
	vec4 texel = get_texture_val(normalize(world_normal), vpos);
	vec3 color = vec4(0.0);
	vec3 tc    = 0.5*(vpos.zxy + 1.0);
	float light_scales[8] = float[](1,1,1,1,1,1,1,1);
	light_scales[0] = texture3D(shadow_tex, tc).r; // diffuse/shadow term
	light_scales[1] = texture3D(ao_tex,     tc).r; // ambient term

	for (int i = 0; i < num_lights; ++i) { // sun_diffuse, galaxy_ambient, dynamic ...
		color += light_scales[i]*add_pt_light_comp(n, epos, i).rgb;
	}
	gl_FragColor = vec4(texel.rgb * clamp(color, 0.0, 1.0), texel.a * gl_Color.a); // use gl_Color alpha directly
}

