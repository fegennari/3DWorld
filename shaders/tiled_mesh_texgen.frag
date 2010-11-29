uniform sampler2D tex0, tex1, tex2, tex3, tex4, tex5, tex6;

const int NTEX = 5; // sand, dirt, ground, rock, snow
varying float weights[NTEX];

vec4 get_weighted_tex(in vec2 tc) {
	vec4 color = vec4(0,0,0,0); // use a loop here?
	if (weights[0] > 0.0) color += texture2D(tex2, tc) * weights[0];
	if (weights[1] > 0.0) color += texture2D(tex3, tc) * weights[1];
	if (weights[2] > 0.0) color += texture2D(tex4, tc) * weights[2];
	if (weights[3] > 0.0) color += texture2D(tex5, tc) * weights[3];
	if (weights[4] > 0.0) color += texture2D(tex6, tc) * weights[4];
	return color;
}

void main()
{
	// process textures
	//vec4 texel0 = get_weighted_tex(gl_TexCoord[0].st);
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);
	vec4 color = vec4((texel0.rgb * texel1.rgb * gl_Color.rgb), (texel0.a * texel1.a * gl_Color.a));
	gl_FragColor = apply_fog(color); // add fog
}
