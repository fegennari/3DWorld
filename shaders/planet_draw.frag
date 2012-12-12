uniform float atmosphere = 1.0;
uniform vec2 light_scale = vec2(1,1);
uniform vec3 sun_pos, ss_pos;
uniform float sun_radius, ss_radius;
uniform sampler2D tex0;
varying vec4 epos;
varying vec3 normal, world_space_pos;


void main()
{
	vec4 texel     = texture2D(tex0, gl_TexCoord[0].st);
	vec3 norm      = normalize(normal); // renormalize
	float lt_atten = calc_light_atten(epos, 0);

	if (sun_radius > 0.0 && ss_radius > 0.0) {lt_atten *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);}

	vec3 light_dir = normalize(gl_LightSource[0].position.xyz - epos.xyz);
	vec3 half_vect = normalize(light_dir - normalize(epos.xyz)); // Eye + L = -eye_space_pos + L
	vec3 ambient   = (light_scale[0] * gl_LightSource[0].ambient.rgb * lt_atten) + (light_scale[1] * gl_LightSource[1].ambient.rgb);
	vec3 diffuse   = gl_LightSource[0].diffuse.rgb * max(dot(norm, light_dir), 0.0) * lt_atten;
	vec3 specular  = gl_FrontLightProduct[0].specular.rgb * pow(max(dot(norm, half_vect), 0.0), gl_FrontMaterial.shininess) * pow(texel.b, 4.0) * lt_atten;
	vec3 color     = (texel.rgb * (ambient + diffuse)) + specular;
	float cloud_val= atmosphere*gen_cloud_alpha(world_space_pos);
	if (cloud_val > 0.0) {color = cloud_val*(ambient + diffuse) + (1.0 - cloud_val)*color;} // no clouds over high mountains?
	gl_FragColor = apply_fog(gl_Color * vec4((color + gl_FrontMaterial.emission.rgb), 1.0));
}
