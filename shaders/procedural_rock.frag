varying vec3 normal;
varying vec4 epos;

void main()
{
	vec3 n     = normalize(normal);
	vec3 color = vec3(0.0);

	for (int i = 0; i < num_lights; ++i) { // sun_diffuse, galaxy_ambient, dynamic ...
		color += add_pt_light_comp(n, epos, i).rgb;
	}
	gl_FragColor = apply_fog_epos(vec4(color, 1.0), epos); // apply standard fog
}

