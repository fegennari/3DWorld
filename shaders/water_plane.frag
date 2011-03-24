varying vec3 normal;
uniform sampler2D tex0;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec3 normal2 = normalize(normal); // renormalize
	vec4 color = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal2, 0);
	if (enable_light1) color += add_light_comp(normal2, 1);
	vec4 frag_color = vec4(texel.rgb * color.rgb, texel.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	//float fog = clamp((gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
	//frag_color.a = mix(1.0, frag_color.a, fog);
	gl_FragColor = apply_fog(frag_color);
}
