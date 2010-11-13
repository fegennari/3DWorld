vec4 add_light_comp(in vec3 normal, in int i) {

	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(vec3(gl_LightSource[i].position));
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	
	// compute the lighting
	vec4 diffuse = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * gl_LightSource[i].ambient;
	return ambient + max(dot(normal, lightDir), 0.0)*diffuse;
}


void main()
{
	//gl_TexCoord[0] = gl_MultiTexCoord0;
	//gl_Position = ftransform();
	gl_Position = gl_Vertex;

	// transform the normal into eye space and normalize
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
	
	// Compute the globalAmbient term
	gl_FrontColor = gl_Color * gl_LightModel.ambient;
	
	if (enable_light0) gl_FrontColor += add_light_comp(normal, 0);
	if (enable_light1) gl_FrontColor += add_light_comp(normal, 1);
} 
