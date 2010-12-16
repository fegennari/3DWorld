uniform sampler2D tex0;
uniform float min_alpha = 0.0;

varying vec3 eye, dlpos, normal; // world space

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	vec3 lit_color = gl_Color.rgb;
	if (enable_dlights) lit_color += add_dlights(dlpos, normal, eye); // dynamic lighting
	vec4 color = vec4(texel.rgb * lit_color.rgb, texel.a * gl_Color.a);
	gl_FragColor = apply_fog(color);
}
