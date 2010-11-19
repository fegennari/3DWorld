vec4 add_light_comp(in vec3 normal, in int i) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	
	// compute the ambient and diffuse lighting
	return gl_FrontLightProduct[i].ambient + max(dot(normal, lightDir), 0.0)*gl_FrontLightProduct[i].diffuse;
}

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space, not normalized
	vec4 color = gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
} 
