uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec4 epos;
varying vec3 normal; // world space

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	vec3 n = normalize(normal); // renormalize
	vec4 color = vec4(0,0,0,1);
	color += add_pt_light_comp(n, epos, 0);
	color += add_pt_light_comp(n, epos, 1);
	color += add_pt_light_comp(n, epos, 2);
	color += add_pt_light_comp(n, epos, 3);
	color += add_pt_light_comp(n, epos, 4);
	color += add_pt_light_comp(n, epos, 5);
	color += add_pt_light_comp(n, epos, 6);
	color += add_pt_light_comp(n, epos, 7);
	color = color * gl_Color.rgba + gl_FrontMaterial.emission;
	gl_FragColor = vec4(texel.rgb * color.rgb, texel.a * gl_Color.a); // use diffuse alpha directly
}
