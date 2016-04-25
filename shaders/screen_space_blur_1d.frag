uniform sampler2D frame_buffer_tex;
uniform int dim_val = 0;

const int BLUR_RADIUS = 4;

in vec2 tc;

void main() {
	vec3 color  = vec3(0.0);
	float tot_w = 0.0;

	for (int v = -BLUR_RADIUS; v <= BLUR_RADIUS; ++v) {
		float weight = exp(-abs(v)); // Gaussian
		color       += weight*textureOffset(frame_buffer_tex, tc, ivec2(v*(1 - dim_val), v*dim_val)).rgb;
		tot_w       += weight;
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
