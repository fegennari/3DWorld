#define TEST_CLIP_T(reg, va, vb, vd, vc) {float t = ((va) - (vb))/(vd); if ((vc) > 0.0) \
                                          {if (t > tmin) tmin = t;} else {if (t < tmax) tmax = t;}}

struct pt_pair {
	vec3 v1, v2;
};

pt_pair clip_line(in vec3 v1, in vec3 v2, in float[6] bounds) {
	// clip to scene bounds
	float tmin = 0.0, tmax = 1.0;
	vec3 dv = v2 - v1;
	TEST_CLIP_T(0x01, bounds[0], v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02, bounds[1], v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, bounds[2], v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08, bounds[3], v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10, bounds[4], v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20, bounds[5], v1.z, dv.z, -dv.z); // +z plane
	pt_pair res;
	
	if (tmin >= tmax) { // clipped away
		res.v1 = v1;
		res.v2 = v1;
	}
	else {
		res.v1 = v1 + dv*tmax;
		res.v2 = v1 + dv*tmin;
	}
	return res;
} 
