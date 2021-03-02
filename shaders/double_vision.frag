uniform sampler2D frame_buffer_tex;
uniform float time      = 0.0;
uniform float intensity = 1.0;
uniform float pos_scale = 1.0;

in vec2 tc;

void main() {
	vec2 offset  = vec2(0.05*intensity*(0.5 - abs(tc.x - 0.5)), 0.0); // horizontal offset, highest at the center of the screen
	fg_FragColor = vec4(0.5*(texture(frame_buffer_tex, (tc - offset)).rgb + texture(frame_buffer_tex, (tc + offset)).rgb), 1.0);
}
