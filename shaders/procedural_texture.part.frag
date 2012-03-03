uniform float noise_scale = 1.0;
uniform float tex_mix_saturate = 1.0;
uniform sampler3D noise_tex;

float procedural_eval(in vec3 pos) {
	float val = 0.0;
	val += 1.00*texture3D(noise_tex, 1.00*pos*noise_scale).r;
	//val += 0.50*texture3D(noise_tex, 1.96*pos*noise_scale).r;
	//val += 0.25*texture3D(noise_tex, 4.03*pos*noise_scale).r;
	return clamp(tex_mix_saturate*(val - 0.5), -0.5, 0.5) + 0.5;
}

