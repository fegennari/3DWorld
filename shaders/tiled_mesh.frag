uniform sampler2D tex0, tex1, texw;
varying float water_depth;

void main()
{
	vec4 texel0  = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1  = texture2D(tex1, gl_TexCoord[1].st);
	float uw_light = 2.0*clamp(10.0*water_depth, 0.0, 1.0)*(texture2D(texw, gl_TexCoord[2].st).g - 0.575);
	vec4 color   = clamp(vec4((texel0.rgb * texel1.rgb * gl_Color.rgb + uw_light), (texel0.a * texel1.a * gl_Color.a)), 0.0, 1.0);
	gl_FragColor = apply_fog(color); // add fog
}
