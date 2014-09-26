uniform vec4 emission = vec4(0,0,0,1);
in vec3 normal;

void main()
{
	vec3 normal2 = normalize(normal); // renormalize
	vec3 color = emission.rgb;
	if (enable_light0) color += add_light_comp0(normal2).rgb;
	if (enable_light1) color += add_light_comp1(normal2).rgb;
	fg_FragColor = vec4(color, gl_Color.a);
}
