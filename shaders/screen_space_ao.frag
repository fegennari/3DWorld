uniform float znear, zfar;
uniform vec2 xy_step;
uniform sampler2D depth_tex;
uniform int NUM_DIRS  = 8;
uniform int NUM_STEPS = 4;

in vec2 tc;

float get_linear_depth(in vec2 pos) {
	float d = texture(depth_tex, pos).r;
	return (2.0 * znear) / (zfar + znear - d * (zfar - znear));
}

void main()
{
	//fg_FragColor = vec4(vec3(get_linear_depth(tc)),1); return;
	float depth0   = get_linear_depth(tc) - 0.0001;
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
			float depth = get_linear_depth(pos);
			if (depth + 0.1 < depth0) {break;} // large depth disconuity, skip this dir

			if (depth < depth0) {
				//if (s == 0) {denom -= 1.0; break;} // if first sample is closer, assume this is the edge of the feature and the back of the face and discard
				weight += (1.0 - float(s)*step_mul);
				break;
			}
		}
	}
	float darken = max(0.0, 1.0*weight/max(denom, 1.0)-0.1);
	fg_FragColor = vec4(0.0, 0.0, 0.0, darken); // darken by weight
}
