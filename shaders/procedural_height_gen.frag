uniform float x0, y0, dx, dy, rx, ry, zscale;

in vec3 vpos;

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


