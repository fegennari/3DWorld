varying vec2 tc;

void main()
{
	float val    = 1.0 - 2.0*abs(tc.s - 0.5);
	val          = min(1.0, 1.1*val); // increase centerline max value width
	float alpha  = pow(val, 3.0); // sharper falloff
	if (alpha <= 0.01) discard; // optional
	fg_FragColor = gl_Color*vec4(1,1,1, alpha);
}
