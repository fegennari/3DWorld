vec4 add_light_comp(in vec3 normal, in int i) {

	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
	
	// compute diffuse term
	vec4 diffuse = max(dot(normal, lightDir), 0.0) * gl_FrontLightProduct[i].diffuse;
	
	// compute the specular term
	float NdotHV  = max(dot(normal, normalize(gl_LightSource[i].halfVector.xyz)), 0.0);
	vec4 specular = gl_FrontLightProduct[i].specular * pow(NdotHV, gl_FrontMaterial.shininess);
	return (gl_FrontLightProduct[i].ambient + specular + diffuse)/gl_LightSource[i].constantAttenuation;
}

