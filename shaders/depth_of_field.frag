uniform sampler2D frame_buffer_tex;
uniform vec2 xy_step;
uniform float focus_depth;
uniform float dof_val;
uniform float dim_val;

const int MAX_BLUR_RADIUS = 8;

in vec2 tc;

void main() {

	vec3 color    = vec3(0.0);
	float tot_w   = 0.0;
	float depth   = get_linear_depth_zval(tc);
	float falloff = clamp(dof_val/max(0.001, abs(depth - focus_depth)), 0.1, 1000.0);

	for (int v = -MAX_BLUR_RADIUS; v <= MAX_BLUR_RADIUS; ++v) {
		float weight = exp(-falloff*abs(v)); // Gaussian
		vec2 pos     = tc + vec2(v*(1.0 - dim_val), v*dim_val)*xy_step;
		color       += weight*texture(frame_buffer_tex, pos).rgb;
		tot_w       += weight;
	}
	color /= tot_w;
	fg_FragColor = vec4(color, 1.0);
}
