uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;

varying vec3 eye, vpos; // world space

const float SMOKE_SCALE = 0.25;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	
	if (eye == vpos) {
		if (color.a == 0.0) discard;
		gl_FragColor = color;
		return;
	}
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
	vec3 pos   = (vpos - off)/scale;
	vec3 dir   = eye - vpos;
	vec3 delta = normalize(dir)*step_delta/scale;
	float nsteps = length(dir)/step_delta;
	int num_steps = 1 + min(1000, int(nsteps)); // round up
	float step_weight = fract(nsteps);
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec3 p = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 tex_val = texture3D(smoke_tex, p.zxy); // rgba = {color.rgb, smoke}
		float smoke = SMOKE_SCALE*tex_val.a*step_weight;
		vec3 rgb_comp = (tex_val.rgb * gl_Fog.color.rgb);
		color = ((color.a == 0.0) ? vec4(rgb_comp, smoke) : mix(color, vec4(rgb_comp, 1.0), smoke));
		pos += delta*step_weight;
		step_weight = 1.0;
	}
	if (color.a == 0.0) discard;
	gl_FragColor = color;
	// update depth if density is large?
}
