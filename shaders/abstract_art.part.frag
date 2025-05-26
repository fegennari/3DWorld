
vec3 colorize(float val) {
	float a = 5*val, b = 7*val, c = 11*val;
	return vec3((a - int(a)), (b - int(b)), (c - int(c)));
}
vec2 conjugate   (vec2 v) {return vec2((v.x*v.x - v.y*v.y), -2.0*v.x*v.y);}
vec2 complex_mult(vec2 v) {return vec2((v.x*v.x - v.y*v.y), (v.x*v.y + v.y*v.x));}

vec4 gen_abstract_art(vec2 tc, vec3 seed) {
	// seed.r selects the fractal, seed.g selects the location, and seed.b selects the zoom level
	int mode = int(3.0*seed.r);
	int cix  = min(4, int(4.0*seed.g));
	vec2 pos = 2.0*tc - vec2(1.0); // [-1.0, 1.0]
	pos.y    = -pos.y; // invert Y
	vec2 c; // center of window

	if (mode == 0) { // madelbrot
		vec2 centers[4] = {vec2(-1.41117, 0.0), vec2(-1.3953, 0.018367), vec2(-1.23102, 0.167902), vec2(0.437802, 0.340398)};
		c = centers[cix];
	}
	else if (mode == 1) { // tricorn
		vec2 centers[4] = {vec2(-1.41975, 0.0), vec2(0.7093, 1.22854), vec2(0.742104, 1.2854), vec2(-1.36524, -0.0388551)};
		c = centers[cix];
	}
	else { // burning ship
		vec2 centers[4] = {vec2(-1.57535, -0.00717238), vec2(0.807586, -1.40662), vec2(-0.676457, -1.10419), vec2(0.810333, -1.40356)};
		c = centers[cix];
	}
	c       += (0.001 + 0.009*seed.b)*pos;
	vec2 z   = vec2(0.0, 0.0);
	uint val = 0;

	for (; val < 200; ++val) {
		if (dot(z, z) > 4.0) break;
		if      (mode == 0) {z = complex_mult    (z)  + c;} // madelbrot
		else if (mode == 1) {z = conjugate       (z)  + c;} // tricorn
		else                {z = complex_mult(abs(z)) + c;} // burning ship
	}
	float val2 = (float(val) - log2(log2(dot(z, z))) + 1.0)/200.0; // from http://www.iquilezles.org/www/articles/mset_smooth/mset_smooth.htm
	//float val2 = val/200.0;
	return vec4(colorize(val2), 1.0);
}

