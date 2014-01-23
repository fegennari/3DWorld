uniform sampler2D tex0, alpha_mask_tex;
uniform float min_alpha = 0.0;
varying vec4 epos;
varying vec3 normal; // eye space

void main()
{
#ifdef ALPHA_MASK_TEX
	if (min_alpha != 0.0) { // this test may or may not help performance
		float alpha_mask = texture2D(alpha_mask_tex, gl_TexCoord[0].st).a;
		if (alpha_mask < min_alpha) discard; // slow
	}
#endif
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	//if (texel.a <= min_alpha) discard; // slow
	vec3 n = (gl_FrontFacing ? normalize(normal) : -normalize(normal)); // two-sided lighting
	vec4 color = vec4(0);
#ifdef SHADOW_ONLY_MODE
	color += add_pt_light_comp(n, epos, 0); // light 0 is the system light
#else
	for (int i = 0; i < 8; ++i) {
		color += add_pt_light_comp(n, epos, i);
	}
	color += gl_FrontMaterial.emission;
#endif
	gl_FragColor = vec4(texel.rgb * clamp(color.rgb, 0.0, 1.0), texel.a * gl_Color.a); // use diffuse alpha directly
}
