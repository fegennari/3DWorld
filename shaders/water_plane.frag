varying vec3 normal;
varying vec4 epos, proj_pos;
uniform sampler2D water_tex, reflection_tex;
uniform vec4 water_color, reflect_color;

void main()
{
	vec4 color = texture2D(water_tex, gl_TexCoord[0].st) * water_color;
	vec3 norm  = normalize(normal); // renormalize

	// calculate lighting
	vec4 lighting = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) lighting += add_light_comp_pos(norm, epos, 0);
	if (enable_light1) lighting += add_light_comp_pos(norm, epos, 1);

	if (reflections) {
		// add some green at shallow view angles
		vec3 epos_n = normalize(epos.xyz);
		color = mix(color, vec4(0.0, 1.0, 0.5, color.a), 0.2*(1.0 - abs(dot(epos_n, norm))));

		// calculate reflections
		float reflect_w  = get_fresnel_reflection(-1.0*epos_n, norm, 1.0, 1.333);
		vec2 ref_tex_st  = 0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5);
		vec4 reflect_tex = vec4(texture2D(reflection_tex, ref_tex_st).rgb, 1.0);
		color = mix(color, reflect_color * reflect_tex, reflect_w);
	}

	// determine final color with fog
	vec4 frag_color = vec4(color.rgb * lighting.rgb, color.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	gl_FragColor = apply_fog(frag_color);
}
