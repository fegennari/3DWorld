uniform sampler2D frame_buffer_tex;
uniform float time      = 0.0;
uniform float intensity = 1.0;

in vec2 tc;

void main() {
	vec3 color  = texture(frame_buffer_tex, tc).rgb;
	vec3 effect = fract(0.008*time*vec3(1.0, 1.1, 1.2) + color.gbr);
	if      (effect.r > max(effect.g, effect.b)) {effect.gb *= 0.5;} // less pastel
	else if (effect.g > max(effect.r, effect.b)) {effect.rb *= 0.5;}
	else if (effect.b > max(effect.g, effect.r)) {effect.gr *= 0.5;}
	effect *= max(min((color.r + color.g + color.b), 1.0), 0.5); // not too bright if background is dark
	fg_FragColor = vec4(mix(color, effect, intensity), 1.0);
}
