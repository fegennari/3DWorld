uniform sampler2D tex0, tex1, tex2, tex3, tex4, tex5, tex6;
uniform float ts2, ts3, ts4, ts5, ts6; // texture scales
uniform float cs2, cs3, cs4, cs5, cs6; // color scales

void main()
{
	vec2 tc = gl_TexCoord[0].st;
	vec4 weights = texture2D(tex0, tc);
	float weights4 = (1.0 - weights.r - weights.g - weights.b - weights.a);
	vec3 texel0  = cs2*weights.r*texture2D(tex2, ts2*tc).rgb +
	               cs3*weights.g*texture2D(tex3, ts3*tc).rgb +
				   cs4*weights.b*texture2D(tex4, ts4*tc).rgb +
				   cs5*weights.a*texture2D(tex5, ts5*tc).rgb +
				   cs6*weights4 *texture2D(tex6, ts6*tc).rgb;
	vec3 texel1  = texture2D(tex1, gl_TexCoord[1].st).rgb; // detail texture
	vec4 color   = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), gl_Color.a);
	gl_FragColor = apply_fog(color); // add fog
}
