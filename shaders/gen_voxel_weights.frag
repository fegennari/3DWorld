uniform vec3 offset, scale;
uniform float start_mag, start_freq, rx, ry, zval;

in vec3 vpos;

vec3 mod289(in vec3 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec4 mod289(in vec4 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec4 permute(in vec4 x) {
	return mod289(((x * 34.0) + 1.0) * x);
}

float simplex(in vec3 v)
{
	vec2 C = vec2(1.0 / 6.0, 1.0 / 3.0);
	vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

	// First corner
	vec3 i = floor(v + dot(v, vec3(C.y)));
	vec3 x0 = v - i + dot(i, vec3(C.x));

	// Other corners
	vec3 g = step(vec3(x0.y, x0.z, x0.x), x0);
	vec3 l = vec3(1.0 - g);
	vec3 i1 = min(g, vec3(l.z, l.x, l.y));
	vec3 i2 = max(g, vec3(l.z, l.x, l.y));

	//   x0 = x0 - 0.0 + 0.0 * C.xxx;
	//   x1 = x0 - i1  + 1.0 * C.xxx;
	//   x2 = x0 - i2  + 2.0 * C.xxx;
	//   x3 = x0 - 1.0 + 3.0 * C.xxx;
	vec3 x1 = x0 - i1 + C.x;
	vec3 x2 = x0 - i2 + C.y; // 2.0*C.x = 1/3 = C.y
	vec3 x3 = x0 - D.y;      // -1.0+3.0*C.x = -0.5 = -D.y

	// Permutations
	i = mod289(i);
	vec4 p = vec4(permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0)));

	// Gradients: 7x7 points over a square, mapped onto an octahedron.
	// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
	float n_ = 0.142857142857; // 1.0/7.0
	vec3 ns = n_ * vec3(D.w, D.y, D.z) - vec3(D.x, D.z, D.x);

	vec4 j = vec4(p - 49.0 * floor(p * ns.z * ns.z));  //  mod(p,7*7)
	vec4 x_ = vec4(floor(j * ns.z));
	vec4 y_ = vec4(floor(j - 7.0 * x_));    // mod(j,N)
	vec4 x = vec4(x_ * ns.x + ns.y);
	vec4 y = vec4(y_ * ns.x + ns.y);
	vec4 h = vec4(1.0 - abs(x) - abs(y));
	vec4 b0 = vec4(x.x, x.y, y.x, y.y);
	vec4 b1 = vec4(x.z, x.w, y.z, y.w);

	// vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
	// vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
	vec4 s0 = vec4(floor(b0) * 2.0 + 1.0);
	vec4 s1 = vec4(floor(b1) * 2.0 + 1.0);
	vec4 sh = -step(h, vec4(0.0));

	vec4 a0 = vec4(b0.x, b0.z, b0.y, b0.w) + vec4(s0.x, s0.z, s0.y, s0.w) * vec4(sh.x, sh.x, sh.y, sh.y);
	vec4 a1 = vec4(b1.x, b1.z, b1.y, b1.w) + vec4(s1.x, s1.z, s1.y, s1.w) * vec4(sh.z, sh.z, sh.w, sh.w);
	vec3 p0 = vec3(a0.x, a0.y, h.x);
	vec3 p1 = vec3(a0.z, a0.w, h.y);
	vec3 p2 = vec3(a1.x, a1.y, h.z);
	vec3 p3 = vec3(a1.z, a1.w, h.w);

	// Normalise gradients
	//vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
	vec4 norm = 1.79284291400159 - 0.85373472095314 * vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3));
	p0 *= norm.x;
	p1 *= norm.y;
	p2 *= norm.z;
	p3 *= norm.w;

	// Mix final noise value
	vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), vec4(0));
	m = m * m;
	return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

void main()
{
	vec3 pos   = offset + scale*vec3(vpos.x, vpos.y, zval);
	float val  = 0.0;
	float mag  = start_mag;
	float freq = start_freq;
	float crx  = rx;
	float cry  = ry;
	const float lacunarity = 1.92;
	const float gain       = 0.5;
	const int NUM_OCTAVES  = 5;

	for (int i = 0; i < NUM_OCTAVES; ++i) {
		float noise = simplex(freq*pos + vec3(crx, cry, crx-cry));
		val  += mag*noise;
		mag  *= gain;
		freq *= lacunarity;
		crx  *= 1.5;
		cry  *= 1.5;
	}
	fg_FragColor = vec4(val, 0.0, 0.0, 1.0); // only red channel is used
}


