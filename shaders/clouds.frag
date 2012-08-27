uniform sampler2D tex0;
uniform float time = 0.0;
uniform float min_alpha = 0.0;

varying vec4 pos, color;

void main()
{
	float alpha = 0.0;
	float freq = 1.0;
	
	for (int i = 0; i < 10; ++i) {
		alpha += texture2D(tex0, 0.001*freq*pos.xy).r/freq; // FIXME: add time
		freq *= 2.0;
	}
	alpha = clamp(4.0*(alpha-1.7), 0.0, 1.0);
	gl_FragColor = gl_FragColor = apply_fog(color*vec4(1,1,1, alpha));
}
