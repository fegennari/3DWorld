void main() {
	if (length(gl_PointCoord - vec2(0.5)) > 0.5) discard; // clip to a circle
	fg_FragColor = gl_Color;
}
