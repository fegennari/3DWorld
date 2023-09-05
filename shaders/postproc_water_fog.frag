uniform sampler2D frame_buffer_tex;
uniform float water_atten;
uniform vec3 uw_atten_max, uw_atten_scale;

in vec2 tc;

void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*atten); // apply scattering and absorption; faster and more tweak-able than using exp()
}

void main() {

	vec3 color   = texture(frame_buffer_tex, tc).rgb;
	float depth  = get_linear_depth_zval(tc);
	atten_color(color, depth*water_atten);
	fg_FragColor = vec4(color, 1.0);
}
