uniform float noise_scale = 1.0;
uniform float tex_mix_saturate = 1.0;
uniform vec3 tex_eval_offset = vec3(0,0,0);
uniform sampler3D noise_tex;

float procedural_eval(in vec3 pos) {
	//return clamp(1.5*abs(sin((pos + 8.0*texture3D(noise_tex, 0.04*pos).r).y)), 0.0, 1.0);
	//return pow(clamp(1.1*sin((pos + 8.0*texture3D(noise_tex, 0.04*pos).r).y), 0.0, 1.0), 3.0);
	float val = 0.0;
	//pos += 8.0*texture3D(noise_tex, 0.004*pos).r; // warp to add veins
	vec3 tc = (pos + tex_eval_offset).zxy*noise_scale;
	val += 1.00*texture3D(noise_tex, 1.00*tc).r;
	//val += 0.50*texture3D(noise_tex, 1.96*tc).r;
	//val += 0.25*texture3D(noise_tex, 4.03*tc).r;
	return clamp(tex_mix_saturate*(val - 0.5), -0.5, 0.5) + 0.5;
}

