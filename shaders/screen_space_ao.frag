uniform int NUM_DIRS   = 8;
uniform int NUM_STEPS  = 4;
uniform float step_len = 4.0; // in pixels

in vec2 tc;

void main() {
	//fg_FragColor = vec4(vec3(get_linear_depth_01(tc)),1); return;
	float depth0   = get_linear_depth_01(tc) - 0.0001;
	//if (depth0 > 0.25) discard; // skip SSAO for high depth (sky/background)
	float dir_mul  = 2.0 * 3.14159 / NUM_DIRS;
	float step_mul = 1.0 / NUM_STEPS;
	float weight   = 0.0;
	float denom    = NUM_DIRS;

	// http://john-chapman-graphics.blogspot.ca/2013/01/ssao-tutorial.html
	//vec3 normal = normalize(cross(dFdy(position), dFdx(position))); 
	
	for (int d = 0; d < NUM_DIRS; d++) {
		vec2 pos    = tc;
		float theta = d*dir_mul;
		vec2 dir    = vec2(sin(theta), cos(theta));
		vec2 step   = step_len * xy_step * dir;

		for (int s = 0; s < NUM_STEPS; s++) {
			pos += step;
			float depth = get_linear_depth_01(pos);
			if (depth + 0.005 < depth0) {break;} // large depth disconuity, skip this dir

			if (depth < depth0) {
				//if (s == 0) {denom -= 1.0; break;} // if first sample is closer, assume this is the edge of the feature and the back of the face and discard
				weight += (1.0 - float(s)*step_mul);
				break;
			}
		}
	}
	float darken = clamp(2.0*(weight/max(denom, 1.0)-0.5), 0.0, 1.0);
#ifdef WRITE_COLOR
    fg_FragColor = vec4(1.0-darken, 1.0-darken, 1.0-darken, 1.0); // write grayscale values
#else
    fg_FragColor = vec4(0.0, 0.0, 0.0, darken); // write as alpha blend of black color
#endif
}
