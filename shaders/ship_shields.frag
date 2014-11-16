uniform sampler2D tex0;
uniform sampler3D noise_tex;
uniform float min_alpha   = 0.0;
uniform float noise_scale = 0.05;
uniform float time        = 0.0;
#define NUM_OCTAVES 6

in vec2 tc;
in vec3 vertex;

void main()
{
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;
	vec3 time_v = 0.04*time*vec3(1.0, 1.2, 1.3);
	float val  = 0.0;
	float freq = 1.0;

	for (int i = 0; i < NUM_OCTAVES; ++i) { // use highly ridged noise
		float v = texture(noise_tex, noise_scale*(freq*vertex + time_v)).r;
		v = 2.0*v - 1.0; // map [0,1] range to [-1,1]
		v = 1.0 - abs(v); // ridged noise
		val  += pow(v, 5.0)/freq;
		freq *= 2.0;
	}
	texel.a *= 0.2 + 1.8*clamp(1.2*(val-0.4), 0.0, 1.0);
	if (texel.a <= min_alpha) discard;
	fg_FragColor = gl_Color * texel;
}
