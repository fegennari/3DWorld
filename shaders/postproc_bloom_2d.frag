uniform sampler2D frame_buffer_tex;
uniform vec2 xy_step;

const int BLOOM_RADIUS = 8;

in vec2 tc;

void main() {
	vec3 base_color  = texture(frame_buffer_tex, tc).rgb;
	vec3 bloom_color = vec3(0.0);
	float falloff    = 0.1;
	float min_val    = 0.2;
	float tot_w      = 0.0;
	float mult       = 2.0/(1.0 - min_val);

	for (int vy = -BLOOM_RADIUS; vy <= BLOOM_RADIUS; ++vy) {
		for (int vx = -BLOOM_RADIUS; vx <= BLOOM_RADIUS; ++vx) {
			float weight = exp(-falloff*length(vec2(vx, vy))); // Gaussian
			vec2 pos     = tc + vec2(vx, vy)*xy_step;
			vec3 color   = weight*texture(frame_buffer_tex, pos).rgb;
			float val    = mult*(max(max(color.r, color.g), color.b) - min_val); // max across all colors
			bloom_color += weight*color*clamp(val, 0.0, 1.0);
			tot_w       += weight;
		}
	}
	bloom_color *= 4.0/tot_w;
	fg_FragColor = vec4(min(vec3(1.0), (base_color + bloom_color)), 1.0);
}
