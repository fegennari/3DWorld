uniform sampler2D frame_buffer_tex;
uniform vec2 xy_step;

const int BLUR_RADIUS = 4;

in vec2 tc;

void main() {
	float color = 0.0; // grayscale
	float tot_w = 0.0;

	for (int y = -BLUR_RADIUS; y <= BLUR_RADIUS; ++y) {
		for (int x = -BLUR_RADIUS; x <= BLUR_RADIUS; ++x) {
			float weight = exp(-sqrt(x*x + y*y)); // Gaussian
			vec2 pos     = tc + vec2(x, y)*xy_step;
			color       += weight*texture(frame_buffer_tex, pos).r; // assume grayscale and use only the red component
			tot_w       += weight;
		}
	}
	color /= tot_w;
	fg_FragColor = vec4(0.0, 0.0, 0.0, 1.0-color); // color=0 => make black, color=1 => make transparent
}
