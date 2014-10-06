uniform float x0, y0, dx, dy, rx, ry, zscale;

in vec3 vpos;

#if 0
float swissTurbulence(float2 p, float seed, int octaves, float lacunarity=2.0, float gain=0.5, float warp=0.15)
{
     float sum  = 0;
     float freq = 1.0;
	 float amp  = 1.0;
     vec2 dsum  = vec2(0,0);

     for (int i = 0; i < octaves; ++i) {
         vec3 n = perlinNoiseDeriv((p + warp * dsum)*freq, seed + i);
         sum   += amp * (1 - abs(n.x));
         dsum  += amp * n.yz * -n.x;
         freq  *= lacunarity;
         amp   *= gain * clamp(sum, 0.0, 1.0);
    }
    return sum;
}
#endif

vec3 mod289(in vec3 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec3 permute(in vec3 x) {
	return mod289(((x * 34.0) + 1.0) * x);
}

float simplex(in vec2 v)
{
	vec4 C = vec4(
		 0.211324865405187,  // (3.0 -  sqrt(3.0)) / 6.0
		 0.366025403784439,  //  0.5 * (sqrt(3.0)  - 1.0)
		-0.577350269189626,	 // -1.0 + 2.0 * C.x
		 0.024390243902439); //  1.0 / 41.0

	// First corner
	vec2 i  = floor(v + dot(v, C.yy));
	vec2 x0 = v - i   + dot(i, C.xx);

	// Other corners
	vec2 i1 = (x0.x > x0.y) ? vec2(1, 0) : vec2(0, 1);
	vec4 x12 = x0.xyxy + C.xxzz;
	x12 = vec4(x12.xy - i1, x12.z, x12.w);

	// Permutations
	i = mod(i, vec2(289)); // Avoid truncation effects in permutation
	vec3 p = permute(permute(i.y + vec3(0, i1.y, 1)) + i.x + vec3(0, i1.x, 1));
	vec3 m = max(vec3(0.5) - vec3(dot(x0, x0), dot(x12.xy, x12.xy), dot(x12.zw, x12.zw)), vec3(0));
	m = m * m;
	m = m * m;

	// Gradients: 41 points uniformly over a line, mapped onto a diamond.
	// The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
	vec3 x = 2.0 * fract(p * C.w) - 1.0;
	vec3 h = abs(x) - 0.5;
	vec3 ox = floor(x + 0.5);
	vec3 a0 = x - ox;

	// Normalise gradients implicitly by scaling m
	// Inlined for speed: m *= taylorInvSqrt(a0*a0 + h*h);
	m *= 1.79284291400159 - 0.85373472095314 * (a0 * a0 + h * h);

	// Compute final noise value at P
	vec3 g;
	g.x  = a0.x  * x0.x   + h.x  * x0.y;
	g.yz = a0.yz * x12.xz + h.yz * x12.yw;
	return 130.0 * dot(m, g);
}

void main()
{
	float x = x0 + dx*vpos.x;
	float y = y0 + dy*vpos.y;
	float zval = 0.0;
	float mag  = 1.0;
	float freq = 1.0;
	float crx  = rx;
	float cry  = ry;
	const float lacunarity = 1.92;
	const float gain       = 0.5;
	const int NUM_OCTAVES  = 9;

	for (int i = 0; i < NUM_OCTAVES; ++i) {
		vec2 pos    = vec2((freq*x + crx), (freq*y + cry));
		float noise = simplex(pos);
#ifdef RIDGED
		noise = 0.45 - abs(noise); // ridged
#endif
#ifdef BILLOWY
		noise = abs(noise) - 0.40; // billowy
#endif
		zval += mag*noise;
		mag  *= gain;
		freq *= lacunarity;
		crx  *= 1.5;
		cry  *= 1.5;
	}
	fg_FragColor = vec4(zscale*zval, 0.0, 0.0, 1.0); // only red channel is used
}


