uniform vec3 light_color = vec3(1.0);

void main() {
	fg_FragColor = gl_Color * vec4(light_color, 1.0);
}
