varying vec3 normal;
varying vec4 epos;

void main()
{
	vec3 n     = normalize(normal);
	vec4 color = gl_FrontMaterial.emission;

	for (int i = 0; i < num_lights; ++i) { // sun_diffuse, galaxy_ambient, dynamic ...
		color += add_pt_light_comp(n, epos, i);
	}
	gl_FragColor = apply_fog(color); // apply standard fog
}

