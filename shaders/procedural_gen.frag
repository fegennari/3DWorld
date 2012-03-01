uniform float half_dxy;
uniform sampler2D tex0;
uniform sampler3D smoke_tex; // for indirect lighting
uniform float min_alpha = 0.0;

uniform vec3 const_indir_color = vec3(0,0,0);

// clipped eye position, clipped vertex position, starting vertex position
varying vec3 eye, vpos, spos, normal; // world space

void main()
{
	vec4 texel  = texture2D(tex0, gl_TexCoord[0].st); // *** FIXME: texturing ***
	float alpha = gl_Color.a;
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
	lit_color += gl_FrontMaterial.diffuse.rgb * const_indir_color; // add constant indir
	
	if (indir_lighting) {
		vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
		vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
		vec3 sp    = clamp((spos  - off)/scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec3 indir_light = texture3D(smoke_tex, sp.zxy).rgb; // add indir light color from texture
		lit_color += gl_FrontMaterial.diffuse.rgb * indir_light; // indirect lighting
	}
	if (direct_lighting) { // directional light sources with no attenuation (Note: could add other lights later)
		vec3 n = normalize(eye_norm);
		if (enable_light0) lit_color += add_light_comp_pos_smap(n, epos, 0).rgb;
		if (enable_light1) lit_color += add_light_comp_pos_smap(n, epos, 1).rgb;
	}
	if (enable_dlights) lit_color += gl_FrontMaterial.diffuse.rgb * add_dlights(vpos, normalize(normal), eye); // dynamic lighting
	vec4 color  = vec4((texel.rgb * lit_color), (texel.a * alpha));
	if (color.a <= min_alpha) discard;
	gl_FragColor = apply_fog(color); // apply standard fog
}
