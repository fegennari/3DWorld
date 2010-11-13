vec4 add_light_comp(in bool shadowed, in vec3 normal, in vec4 eye_space_pos, in int i) {

	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(vec3(gl_LightSource[i].position));
	
	// compute the ambient and diffuse terms
	vec4 diffuse = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * gl_LightSource[i].ambient;
	vec4 color = ambient;
	
	if (!shadowed) { // calculate specular lighting and normal effects
		vec3 dir_to_camera = normalize(vec3(0,0,0) - eye_space_pos.xyz);
		vec3 dir_to_light = normalize(gl_LightSource[i].position.xyz - eye_space_pos.xyz);
		float dp1 = dot(normal, dir_to_camera);
		float dp2 = dot(normal, dir_to_light);
		
		if ((dp1 < 0.0) != (dp2 < 0.0)) { // looking at unlit side
			normal *= -0.5; // reverse and halve
			
			if (i == 0) {
				float dp3 = dot(dir_to_camera, dir_to_light);

				if (dp3 < -0.95) { // leaf between light source and eye
					float val = -20.0*(dp3 + 0.95);
					normal = normal*(1.0 - val) + dir_to_light*((dp2 < 0.0) ? 1.0 : -1.0)*val; // max light
				}
			}
		}
		
		// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
		float NdotL = dot(normal, lightDir);
		
		// compute the specular and diffuse terms if not shadowed and NdotL is larger than zero
		if (NdotL > 0.0) {
			float NdotHV = max(dot(normal, normalize(gl_LightSource[i].halfVector.xyz)), 0.0);
			color += gl_FrontMaterial.specular * gl_LightSource[i].specular * pow(NdotHV, gl_FrontMaterial.shininess);
		}
	}
	return color + max(dot(normal, lightDir), 0.0)*diffuse;
}


void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = ftransform();
	bool shadowed = (sqrt(dot(gl_Normal, gl_Normal)) < 0.4);

	// transform the normal into eye space, but don't normalize because it may be scaled for shadows, and we know there are no glScale() calls
	vec3 normal = gl_NormalMatrix * gl_Normal;
	
	vec4 eye_space_pos = gl_ModelViewMatrix * gl_Vertex;
	if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse
	
	// Compute the globalAmbient term
	gl_FrontColor = gl_Color * gl_LightModel.ambient;
	
	if (enable_light0) gl_FrontColor += add_light_comp(shadowed, normal, eye_space_pos, 0);
	if (enable_light1) gl_FrontColor += add_light_comp(shadowed, normal, eye_space_pos, 1);
} 
