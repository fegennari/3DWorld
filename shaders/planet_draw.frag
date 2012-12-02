uniform vec2 light_scale = vec2(1,1);
uniform sampler2D tex0;
varying vec4 epos;
varying vec3 normal; // world space


vec4 add_light_planet(in vec3 normal, in vec4 epos, in vec4 texel, in int i)
{
	// normalize the light's direction in eye space
	vec3 light_dir = normalize(gl_LightSource[i].position.xyz - epos.xyz);
	vec4 diffuse = gl_Color * texel * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * texel * gl_LightSource[i].ambient;
	vec3 half_vect = normalize(light_dir - normalize(epos.xyz)); // Eye + L = -eye_space_pos + L
	vec4 specular  = gl_FrontLightProduct[i].specular * pow(max(dot(normal, half_vect), 0.0), gl_FrontMaterial.shininess)*pow(texel.b, 4.0);
	return (ambient + (max(dot(normal, light_dir), 0.0)*diffuse + specular)) * calc_light_atten(epos, i);
}

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec3 norm  = normalize(normal); // renormalize
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += light_scale[0]*add_light_planet(norm, epos, texel, 0).rgb;
	color.rgb += light_scale[1]*add_light_planet(norm, epos, texel, 1).rgb;
	gl_FragColor = apply_fog(color);
}
