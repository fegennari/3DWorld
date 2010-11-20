varying float fogFactor;
uniform sampler2D tex0, tex1;

void main()
{
	// process tex0 and tex1
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);
	vec4 color = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), (texel0.a * texel1.a * gl_Color.a));
	
	// add fog
	gl_FragColor = mix(gl_Fog.color, color, fogFactor);
}
