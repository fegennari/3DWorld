uniform sampler2D frame_buffer_tex;
uniform vec4 edge_color = vec4(0.0);

in vec2 tc;

void main() {
	vec2 edge_dist    = min(tc, (vec2(1.0)-tc));
	float edge_weight = (1.0 - min(1.0, 200.0*edge_dist.x*edge_dist.y))*edge_color.a;
	vec3 fb_color     = texture(frame_buffer_tex, tc).rgb;
	fg_FragColor      = vec4(mix(fb_color, edge_color.rgb, edge_weight), 1.0);
}
