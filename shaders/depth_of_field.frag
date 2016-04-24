uniform sampler2D frame_buffer_tex;
uniform float focus_depth;
uniform float dof_val;
uniform int dim_val = 0;

const int MAX_BLUR_RADIUS = 8;

in vec2 tc;

void main() {

	vec3 color    = vec3(0.0);
	float tot_w   = 0.0;
	float depth   = get_linear_depth_zval(tc);
	float falloff = clamp(dof_val/max(0.001, abs(depth - focus_depth)), 0.1, 1000.0);

	for (int v = -MAX_BLUR_RADIUS; v <= MAX_BLUR_RADIUS; ++v) {
		float weight = exp(-falloff*abs(v)); // Gaussian - Note: could use a lookup table, but doesn't make much difference
		color       += weight*textureOffset(frame_buffer_tex, tc, ivec2(v*(1 - dim_val), v*dim_val)).rgb;
		tot_w       += weight;
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
