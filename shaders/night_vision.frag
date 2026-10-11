uniform float xscale=0.0, yscale=0.0;
uniform sampler2D frame_buffer_tex;

in vec2 tc;

void main() {
	vec3 color      = texture(frame_buffer_tex, vec2(tc.x*xscale, tc.y*yscale)).rgb;
	float luminance = (color.r + color.g + color.b)/3.0;
	float gscale    = mix(3.0*luminance, luminance, luminance);
	fg_FragColor    = vec4(0.0, gscale, 0.0, 1.0); // all green
}
