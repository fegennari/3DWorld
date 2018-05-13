uniform sampler2D tex0;
uniform float min_alpha = 0.0;

in vec2 tc;

void main() {
	vec4 texel   = texture(tex0, tc);
	float radius = length(2.0*(tc - vec2(0.5)));
	texel.a     *= clamp(4.0*(1.0 - radius), 0.0, 1.0);
	if (texel.a <= min_alpha) discard;
	fg_FragColor = gl_Color*texel;
}
