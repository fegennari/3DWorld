varying vec4 epos;

void main()
{
	float alpha = 1.0/(1.0 + 50.0*length(epos.xyz)); // attenuate when far from the camera
	if (alpha < 0.04) {discard;}
	vec3 normal = normalize(-epos.xyz); // facing the camera
	vec4 color  = gl_FrontMaterial.emission;

	for (int i = 0; i < 2; ++i) { // sun_diffuse, galaxy_ambient
		color += add_pt_light_comp(normal, epos, i);
	}
	gl_FragColor = vec4(color.rgb, alpha * gl_Color.a); // use diffuse alpha directly;
}

