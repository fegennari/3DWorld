uniform sampler2D tex0;
uniform float min_alpha = 0.0;
uniform vec4  emission  = vec4(0,0,0,1);

varying vec3 dlpos, dl_normal; // world space
varying vec3 normal;
varying vec2 tc;

void main()
{
	vec4 texel = texture2D(tex0, tc);
	if (texel.a <= min_alpha) discard;
	
	vec4 epos    = gl_ModelViewMatrix * vec4(dlpos, 1.0);
	vec3 normal2 = (no_normalize ? normal : normalize(normal)); // renormalize
	vec3 color   = emission.rgb;
	if (enable_dlights) color += add_dlights(dlpos, normalize(dl_normal), vec3(1.0)).rgb; // dynamic lighting
	if (enable_light0 ) color += add_light_comp_pos_smap_light0(normal2, epos).rgb;
	if (enable_light1 ) color += add_light_comp_pos_smap_light1(normal2, epos).rgb;
	vec4 frag_color = vec4(texel.rgb * color, texel.a * gl_Color.a); // use gl_Color alpha directly
#ifndef NO_FOG
	frag_color = apply_fog_epos(frag_color, epos);
#endif
	gl_FragColor = frag_color;
}
