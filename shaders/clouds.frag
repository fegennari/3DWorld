uniform sampler2D tex0;
uniform vec2 dxy;
uniform float min_alpha = 0.0;

varying vec4 pos, color;

void main()
{
	float alpha = 0.0;
	
	for (int n = 0; n < 2; ++n) {
		float freq = 1.0;
		float freq2 = 0.005 * (1.0 + n);

		for (int i = 0; i < 8; ++i) {
			alpha += texture2D(tex0, 0.0015*(freq*pos.xy + freq2*dxy)).r/freq;
			freq  *= 2.0;
			freq2 *= 2.25;
		}
	}
	alpha = clamp(5.0*(0.5*alpha-1.68), 0.0, 1.0);
	gl_FragColor = gl_FragColor = apply_fog(color*vec4(1,1,1, alpha));
}
