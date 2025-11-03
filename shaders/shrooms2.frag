uniform sampler2D frame_buffer_tex;
uniform float time         = 0.0;
uniform float intensity    = 1.0;
uniform float aspect_ratio = 1.0;

in vec2 tc;

float rand_val(float val) {return fract(sin(1.3 * val)*2.7);}

void main() {
	float tval  = 0.001*time;
	float val   = 1.0 + tval;
	vec4 effect = vec4(0.0);

	for (int i = 0; i < 40; ++i) {
		vec2 center  = vec2(rand_val(val), rand_val(val+1.0));
		float dist   = (0.25 + 0.75*rand_val(val+3.0 + 2.3*tval))*length((tc - center)*vec2(aspect_ratio, 1.0));
		float weight = intensity*max(0.0, 40.0*(0.025 - dist));
		vec3 color   = vec3(rand_val(val+2.0), rand_val(val+3.0), rand_val(val+4.0));
		if      (color.r > max(color.g, color.b)) {color.gb *= 0.5;} // less pastel
		else if (color.g > max(color.r, color.b)) {color.rb *= 0.5;}
		else if (color.b > max(color.g, color.r)) {color.gr *= 0.5;}
		else {color.rgb = vec3(1.0);}
		effect += weight*vec4(color, 1.0);
		val    += 5.0 + 0.02*tval;
	}
	float vmax = max(max(0.01, effect.r), max(effect.g, effect.b));
	fg_FragColor = vec4(effect.rgb/vmax, effect.a);
}
