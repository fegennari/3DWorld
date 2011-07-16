uniform sampler2D tex0;
uniform float min_alpha = 0.0;

varying vec4 vertex; // world space
varying vec3 lcolor;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	float cmult = 1.0;
	if (use_shadow_map) cmult = get_shadow_map_weight(vertex, 0);
	vec3 color  = gl_Color.rgb + cmult*lcolor;
	gl_FragColor = apply_fog(vec4(texel.rgb * color, texel.a * gl_Color.a));
}
