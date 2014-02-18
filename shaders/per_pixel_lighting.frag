varying vec3 normal;

void main()
{
	vec3 normal2 = normalize(normal); // renormalize
	vec4 color = gl_FrontMaterial.emission;
	if (enable_light0) color += add_light_comp0(normal2);
	if (enable_light1) color += add_light_comp1(normal2);
	gl_FragColor = color;
}
