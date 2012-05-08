uniform sampler2D tex0, tex1, tex2, tex3, tex4, tex5, tex6;
uniform float ts2, ts3, ts4, ts5, ts6; // texture scales

void main()
{
	vec2 tc = gl_TexCoord[0].st;
	vec4 weights = texture2D(tex0, tc);
	float weights4 = (1.0 - weights.r - weights.g - weights.b - weights.a);
	vec4 texel0  = weights.r*texture2D(tex2, ts2*tc) + weights.g*texture2D(tex3, ts3*tc) + weights.b*texture2D(tex4, ts4*tc) + weights.a*texture2D(tex5, ts5*tc) + weights4*texture2D(tex6, ts6*tc);
	vec4 texel1  = texture2D(tex1, gl_TexCoord[1].st); // detail texture
	vec4 color   = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), (texel0.a * texel1.a * gl_Color.a));
	gl_FragColor = apply_fog_quadratic(color); // add fog
}
