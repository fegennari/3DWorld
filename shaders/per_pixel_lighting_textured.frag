uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space
varying vec3 dl_normal;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	
	vec3 normal2 = (no_normalize ? normal : normalize(normal)); // renormalize
	vec4 color = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_dlights) color.rgb += add_dlights(dlpos, normalize(dl_normal), eye, vec3(1,1,1)).rgb; // dynamic lighting
	if (enable_light0 ) color += add_light_comp_pos_smap(normal2, epos, 0);
	if (enable_light1 ) color += add_light_comp_pos_smap(normal2, epos, 1);
	vec4 frag_color = vec4(texel.rgb * color.rgb, texel.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	gl_FragColor = apply_fog(frag_color);
}
