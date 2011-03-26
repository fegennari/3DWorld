varying vec3 normal;
varying vec4 epos;
uniform sampler2D tex0;
uniform vec4 water_color, reflect_color;

void main()
{
	vec4 color   = texture2D(tex0, gl_TexCoord[0].st) * water_color;
	vec3 normal2 = normalize(normal); // renormalize

	// calculate lighting
	vec4 lighting = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) lighting += add_light_comp_pos(normal2, epos, 0);
	if (enable_light1) lighting += add_light_comp_pos(normal2, epos, 1);

	if (reflections) { // calculate reflections
		float reflect_w = get_fresnel_reflection(-1.0*normalize(epos.xyz), normal2, 1.0, 1.333);
		color = mix(color, reflect_color, reflect_w);
	}

	// determine final color with fog
	vec4 frag_color = vec4(color.rgb * lighting.rgb, color.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	//float fog = clamp((gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
	//frag_color.a = mix(1.0, frag_color.a, fog);
	gl_FragColor = apply_fog(frag_color);
}
