uniform sampler2D tex0, tex1;
varying vec3 eye, dlpos, normal; // world space

void main()
{
	vec4 texel0 = texture2D(tex0, gl_TexCoord[0].st);
	vec4 texel1 = texture2D(tex1, gl_TexCoord[1].st);
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
	if (enable_dlights) lit_color += add_dlights(dlpos, normalize(normal), eye); // dynamic lighting
	vec4 color = vec4((texel0.rgb * texel1.rgb * lit_color), (texel0.a * texel1.a * gl_Color.a));
	gl_FragColor = apply_fog(color); // add fog
}
