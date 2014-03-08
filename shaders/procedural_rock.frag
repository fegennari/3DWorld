varying vec3 normal;
varying vec4 epos;

void main()
{
	vec3 n     = normalize(normal);
	vec4 color = vec4(0,0,0,1);

	for (int i = 0; i < num_lights; ++i) { // sun_diffuse, galaxy_ambient, dynamic ...
		color.rgb += add_pt_light_comp(n, epos, i).rgb;
	}
	gl_FragColor = apply_fog_epos(color, epos); // apply standard fog
}

