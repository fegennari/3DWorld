uniform float x0, y0, dx, dy;

varying vec3 vpos;

void main()
{
	float x = 0.1*x0 + 100.0*dx*vpos.x;
	float y = 0.1*y0 + 100.0*dy*vpos.y;
	float height = abs(sin(x)*sin(y));
	fg_FragColor = vec4(height, 0.0, 0.0, 1.0);
}
