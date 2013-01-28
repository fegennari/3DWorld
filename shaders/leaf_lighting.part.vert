uniform int num_dlights = 0;
uniform float normal_scale = 1.0;
uniform vec4 color_scale = vec4(1,1,1,1);

vec4 add_leaf_light_comp(in bool shadowed, in vec3 normal, in vec4 eye_space_pos, in int i)
{
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
	
	// compute the ambient and diffuse terms
	vec4 diffuse = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * gl_LightSource[i].ambient;
	
	// calculate specular lighting and normal effects
	vec3 dir_to_camera = normalize(vec3(0.0,0.0,0.0) - eye_space_pos.xyz);
	vec3 dir_to_light  = normalize(gl_LightSource[i].position.xyz - eye_space_pos.xyz);
	float dp1 = dot(normal, dir_to_camera);
	float dp2 = dot(normal, dir_to_light);
	
	if (!shadowed && (dp1 < 0.0) != (dp2 < 0.0)) { // looking at unlit side
		normal *= -0.5; // reverse and halve
		float dp3 = dot(dir_to_camera, dir_to_light);

		if (i == 0 && dp3 < -0.95) { // leaf between light source and eye
			float val = -20.0*(dp3 + 0.95);
			normal = normal*(1.0 - val) + dir_to_light*((dp2 < 0.0) ? val : -val); // max light
		}
	}
	
	// compute the specular and diffuse terms if not shadowed
	float NdotHV  = max(dot(normal, normalize(gl_LightSource[i].halfVector.xyz)), 0.0);
	vec4 specular = (shadowed ? 0.0 : 1.0) * gl_FrontMaterial.specular * gl_LightSource[i].specular * pow(NdotHV, gl_FrontMaterial.shininess);
	return (ambient + specular + max(dot(normal, lightDir), 0.0)*diffuse);
}

void calc_leaf_lighting()
{
	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = gl_NormalMatrix * gl_Normal * normal_scale;
	
	vec4 eye_space_pos = gl_ModelViewMatrix * gl_Vertex;
	if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse (could use faceforward())
	
	// Compute the globalAmbient term
	bool shadowed = (sqrt(dot(gl_Normal, gl_Normal)) < 0.4);
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 0);
	if (enable_light1) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 1);
#ifndef NO_LEAF_DLIGHTS
	vec3 n = normalize(normal);

	for (int i = 2; i < num_dlights+2; ++i) { // add N dynamic point light sources
		color += add_pt_light_comp(n, eye_space_pos, i);
	}
#endif
	gl_FrontColor   = clamp(color*color_scale, 0.0, 1.0);
	gl_FogFragCoord = length(eye_space_pos.xyz);
}