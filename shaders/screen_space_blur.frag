uniform sampler2D frame_buffer_tex;
uniform vec2 xy_step;

const int BLUR_RADIUS = 4;

in vec2 tc;

void main() {
	vec3 color  = vec3(0.0);
	float tot_w = 0.0;

	for (int y = -BLUR_RADIUS; y <= BLUR_RADIUS; ++y) {
		for (int x = -BLUR_RADIUS; x <= BLUR_RADIUS; ++x) {
			float weight = exp(-sqrt(x*x + y*y)); // Gaussian
			vec2 pos     = tc + vec2(x, y)*xy_step;
			color       += weight*texture(frame_buffer_tex, pos).rgb;
			tot_w       += weight;
		}
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
