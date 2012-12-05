uniform sampler2D tex0;
varying vec4 epos;
varying vec3 normal; // note: unused

vec4 add_light_rings(in vec3 n, in vec4 eye_pos, in int i)
{
	vec3 light_dir = normalize(gl_LightSource[i].position.xyz - eye_pos.xyz); // normalize the light's direction in eye space
	vec4 diffuse   = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient   = gl_Color * gl_LightSource[i].ambient;
	vec4 specular  = get_light_specular(n, light_dir, eye_pos.xyz, i);
	return (ambient + (abs(dot(n, light_dir))*diffuse + specular)) * calc_light_atten(eye_pos, i);
}

void main()
{
	vec2 tc = 16*gl_TexCoord[0].st;
	vec3 norm2 = normalize(normal + vec3(texture2D(tex0, tc).r-0.5, texture2D(tex0, tc+vec2(0.4,0.7)).r-0.5, texture2D(tex0, tc+vec2(0.3,0.8)).r-0.5));
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += add_light_rings(norm2, epos, 0).rgb;
	color.rgb += add_light_rings(norm2, epos, 1).rgb;
	color.a   *= texture2D(tex0, 25*gl_TexCoord[0].st).r;
	gl_FragColor = apply_fog(color);
}
