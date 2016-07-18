uniform sampler2D frame_buffer_tex;
uniform vec2 xy_step;
uniform int dim_val = 0;

const int BLUR_RADIUS = 4;

in vec2 tc;

void main() {
	vec3 color  = vec3(0.0);
	float tot_w = 0.0;

	for (int v = -BLUR_RADIUS; v <= BLUR_RADIUS; ++v) {
		float weight = exp(-abs(v)); // Gaussian
		vec2 pos     = tc + vec2(v*(1.0 - dim_val), v*dim_val)*xy_step;
		color       += weight*texture(frame_buffer_tex, pos).rgb;
		tot_w       += weight;
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
