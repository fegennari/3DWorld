uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform sampler2D tex0;
uniform float min_alpha = 0.0;

varying vec3 eye, dlpos, normal; // world space

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	if (texel.a <= min_alpha) discard;
	vec3 lit_color = gl_Color.rgb;

	if (enable_dlights) {
		vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
		vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
		vec3 dlp   = clamp((dlpos - off)/scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		lit_color += add_dlights(dlp, off, scale, normal, dlpos, eye, x_scene_size); // dynamic lighting
	}
	vec4 color = vec4(texel.rgb * lit_color.rgb, texel.a * gl_Color.a);
	gl_FragColor = apply_fog(color);
}
