uniform sampler2D cloud_noise_tex;
uniform vec2 dxy;

float gen_alpha_val(in vec2 pos)
{
	float val = 0.0;
	
	for (int n = 0; n < 2; ++n) {
		float freq  = 1.0;
		float freq2 = 0.005 * (1.0 + n);

		for (int i = 0; i < 8; ++i) {
			val   += texture2D(cloud_noise_tex, 0.0015*(freq*pos.xy + freq2*dxy)).r/freq;
			freq  *= 2.0;
			freq2 *= 2.25;
		}
	}
	return clamp(7.0*(0.5*val-1.68), 0.0, 1.0);
}
