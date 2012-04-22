uniform sampler2D tex0;
uniform float min_alpha = 0.0;
varying vec4 epos;
varying vec3 normal; // world space

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	vec3 n = (gl_FrontFacing ? normalize(normal) : -normalize(normal)); // two-sided lighting
	vec4 color = gl_FrontMaterial.emission;

	for (int i = 0; i < 8; ++i) {
		color += add_pt_light_comp(n, epos, i);
	}
	color = vec4(texel.rgb * clamp(color.rgb, 0.0, 1.0), texel.a * gl_Color.a); // use diffuse alpha directly
#ifndef NO_FOG
	color = apply_fog(color); // apply standard fog
#endif
	gl_FragColor = color;
}
