void main() {
	vec2 pos = gl_PointCoord - vec2(0.5);
	if (dot(pos, pos) > 0.25) discard; // clip to a circle
	fg_FragColor = gl_Color;
}
