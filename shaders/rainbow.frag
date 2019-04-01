uniform float alpha_scale = 1.0;

in vec2 tc;

void main() {
	// Note: outer radius is tan(42)*dist, inner radius is tan(40)*dist
	// Note: could add double rainbow at 51 degrees
	float radius = 2.0*length(tc.st - vec2(0.5));
	if (radius > 1.0) discard;
	const float thick = 0.07; // approx (tan(42) - tan(40))/tan(42)
	float val = 8.0*(radius - (1.0 - thick))/thick; // 0-8
	vec4 color = vec4(0);
	if      (val > 6) {color = mix(vec4(1,0,0,1),  vec4(1,0,0,0),  0.5*(val - 6));} // red => transparent red (2x width)
	else if (val > 5) {color = mix(vec4(1,.5,0,1), vec4(1,0,0,1),  (val - 5));} // orange => red
	else if (val > 4) {color = mix(vec4(1,1,0,1),  vec4(1,.5,0,1), (val - 4));} // yellow => orange
	else if (val > 3) {color = mix(vec4(0,1,0,1),  vec4(1,1,0,1),  (val - 3));} // green => yellow
	else if (val > 2) {color = mix(vec4(0,0,1,1),  vec4(0,1,0,1),  (val - 2));} // blue => green
	else if (val > 1) {color = mix(vec4(0,0,.5,1), vec4(0,0,1,1),  (val - 1));} // indigo => blue
	else if (val > 0) {color = mix(vec4(.5,0,1,0), vec4(0,0,.5,1), (val - 0));} // transparent violet => indigo
	else              {color = mix(vec4(1,1,1,0),  vec4(.5,0,1,0), max(0, (val + 1.0)));}
	color += vec4(0.1,0.1,0.1,0.0); // mix in a bit of white
	color.a *= 0.4; // set transparency
	color.a += 0.2*(min(1.0, (1.0 - radius)/thick)); // a bit of white in the center
	color.a *= alpha_scale;
	float depth = get_linear_depth_01(gl_FragCoord.xy*xy_step);
	color.a *= (depth - 0.1)/0.9; // multiply by depth to integrate drolet particle volume over distance
	fg_FragColor = color;
}
