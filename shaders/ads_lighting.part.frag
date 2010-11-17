vec4 add_light_comp(in vec3 normal, int int i) {

	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(vec3(gl_LightSource[i].position));
	
	// compute the ambient and diffuse terms
	vec4 diffuse = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * gl_LightSource[i].ambient;
	
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	
	// compute the specular and diffuse terms if NdotL is larger than zero
	float NdotHV  = max(dot(normal, normalize(gl_LightSource[i].halfVector.xyz)), 0.0);
	vec4 specular = ((NdotL > 0.0) ? 1.0 : 0.0) * gl_FrontMaterial.specular * gl_LightSource[i].specular * pow(NdotHV, gl_FrontMaterial.shininess);
	return (ambient + specular + max(dot(normal, lightDir), 0.0)*diffuse);
}

