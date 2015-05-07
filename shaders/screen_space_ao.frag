uniform vec2 xy_step;
uniform sampler2D depth_tex;
uniform int NUM_DIRS  = 8;
uniform int NUM_STEPS = 4;

in vec2 tc;

void main()
{
	float depth0   = texture(depth_tex, tc).r - 0.0001;
	float dir_mul  = 2.0 * 3.14159 / NUM_DIRS;
	float step_mul = 1.0 / NUM_STEPS;
	float weight   = 0.0;
	float denom    = NUM_DIRS;
	
	for (int d = 0; d < NUM_DIRS; d++) {
		vec2 pos    = tc;
		float theta = d*dir_mul;
		vec2 dir    = vec2(sin(theta), cos(theta));
		vec2 step   = xy_step * dir;

		for (int s = 0; s < NUM_STEPS; s++) {
			pos += step;
			float depth = texture(depth_tex, pos).r;

			if (depth < depth0) {
				//if (s == 0) {denom -= 1.0; break;} // if first sample is closer, assume this is the edge of the feature and the back of the face and discard
				weight += (1.0 - float(s)*step_mul);
				break;
			}
		}
	}
	gl_FragColor = vec4(0.0, 0.0, 0.0, 0.8*weight/denom); // darken by weight
}
