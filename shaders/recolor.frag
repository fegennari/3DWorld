uniform sampler2D frame_buffer_tex;
//uniform sampler2D noise_tex;
//uniform float time      = 0.0;
//uniform float intensity = 1.0;

in vec2 tc;

void main() {
	vec3 c = texture(frame_buffer_tex, tc).rgb;
	//c.r = c.g = c.b = (c.r + c.b + c.g)/3.0; // grayscale
	//c = vec3(1.0) - c; // inverted
	//float temp = c.r; c.r = c.b; c.b = temp; // swap red and blue
	//float gamma = 2.2; float exposure = 1.0; c = clamp(exposure*c, 0.0, 1.0); c = pow(c, vec3(1.0 / gamma)); // gamma correction
	c = abs(sin(8.0*c)); // psycho
	fg_FragColor = vec4(c, 1.0);
}
