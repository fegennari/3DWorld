uniform sampler2D tex0, tex1;
uniform float fog_scale = 0.0;

void main()
{
	// process tex0 and tex1
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);
	vec4 color = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), (texel0.a * texel1.a * gl_Color.a));
	
	// add fog
	float fog = clamp((gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
	gl_FragColor = mix(gl_Fog.color, color, mix(1.0, fog, fog_scale));
}
