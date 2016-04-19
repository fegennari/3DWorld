uniform sampler2D frame_buffer_tex;

const int BLUR_RADIUS = 4;

in vec2 tc;

void main() {
	vec3 color  = vec3(0.0);
	float tot_w = 0.0;

	for (int y = -BLUR_RADIUS; y <= BLUR_RADIUS; ++y) {
		for (int x = -BLUR_RADIUS; x <= BLUR_RADIUS; ++x) {
			float weight = exp(-sqrt(x*x + y*y)); // Gaussian
			color       += weight*textureOffset(frame_buffer_tex, tc, ivec2(x, y)).rgb;
			tot_w       += weight;
		}
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
