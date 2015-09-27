uniform sampler2D frame_buffer_tex;
uniform float intensity = 1.0;
uniform vec2 xy_step;

const int BLOOM_RADIUS = 8;

in vec2 tc;

void main() {
	vec3 base_color  = texture(frame_buffer_tex, tc).rgb;
	vec3 bloom_color = vec3(0.0);
	float falloff    = 0.5;
	float tot_w      = 0.0;

	for (int y = -BLOOM_RADIUS; y <= BLOOM_RADIUS; ++y) {
		for (int x = -BLOOM_RADIUS; x <= BLOOM_RADIUS; ++x) {
			float weight = exp(-falloff*length(vec2(x, y))); // Gaussian
			vec2 pos     = tc + vec2(x, y)*xy_step;
			vec3 color   = texture(frame_buffer_tex, pos).rgb;
			bloom_color += weight*color*clamp(2.0*(color.b - 0.5), 0.0, 1.0);
			tot_w       += weight;
		}
	}
	bloom_color *= 0.8/tot_w;
	fg_FragColor = vec4(min(vec3(1.0), (base_color + bloom_color)), 1.0);
}
