varying vec3 normal;
varying vec4 epos;
uniform sampler2D tex0;
uniform float water_color[4], reflect_color[4];

vec4 unpack_color(in float c[4]) {
	return vec4(c[0], c[1], c[2], c[3]);
}

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st) * unpack_color(water_color);
	vec3 normal2 = normalize(normal); // renormalize

	// calculate lighting
	vec4 color = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos(normal2, epos, 0);
	if (enable_light1) color += add_light_comp_pos(normal2, epos, 1);

	if (reflections) { // calculate reflections
		float reflect_w = get_fresnel_reflection(-1.0*normalize(epos.xyz), normal2, 1.0, 1.333);
		texel = mix(texel, unpack_color(reflect_color), reflect_w);
	}

	// determine final color with fog
	vec4 frag_color = vec4(texel.rgb * color.rgb, texel.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	//float fog = clamp((gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
	//frag_color.a = mix(1.0, frag_color.a, fog);
	gl_FragColor = apply_fog(frag_color);
}
