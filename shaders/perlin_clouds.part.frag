uniform sampler2D cloud_noise_tex;
uniform vec2 dxy;
uniform float cloud_scale;

float gen_cloud_alpha(in vec2 pos)
{
	float alpha = 0.0;
	
	for (int n = 0; n < 2; ++n) {
		float freq  = 1.0;
		float freq2 = 0.005 * (1.0 + n);

		for (int i = 0; i < NUM_OCTAVES; ++i) {
			alpha += texture2D(cloud_noise_tex, 0.0015*(freq*pos.xy + freq2*dxy)).r/freq;
			freq  *= 2.0;
			freq2 *= 2.25;
		}
	}
	return clamp(7.0*(cloud_scale*alpha-1.68), 0.0, 1.0);
}
