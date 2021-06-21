uniform float depth_value;
uniform float min_alpha = 0.0;
uniform sampler2D tex0;

in vec2 tc;

void main() {
#ifndef NO_ALPHA_TEST
	if (texture(tex0, tc).a <= min_alpha) discard;
#endif
	gl_FragDepth = depth_value;
}
