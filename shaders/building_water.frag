uniform sampler2D reflection_tex;

in vec4 proj_pos;

void main() {
	vec2 uv = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	float alpha = 0.5; // TODO: fresnel
	fg_FragColor = vec4(texture(reflection_tex, uv).rgb, alpha);
}
