void main()
{
	vec3 normal = normalized(gl_NormalMatrix * gl_Normal);
	vec4 colorgl_Color * gl_LightModel.ambient; // global ambient
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FragColor = color;
}
