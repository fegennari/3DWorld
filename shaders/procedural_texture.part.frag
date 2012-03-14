uniform float noise_scale = 1.0;
uniform float tex_mix_saturate = 1.0;
uniform vec3 tex_eval_offset = vec3(0,0,0);
uniform sampler3D noise_tex;

float procedural_eval(in vec3 pos) {
	float val = 0.0;
	vec3 tc = (pos.zxy + tex_eval_offset)*noise_scale;
	val += 1.00*texture3D(noise_tex, 1.00*tc).r;
	//val += 0.50*texture3D(noise_tex, 1.96*tc).r;
	//val += 0.25*texture3D(noise_tex, 4.03*tc).r;
	return clamp(tex_mix_saturate*(val - 0.5), -0.5, 0.5) + 0.5;
}

