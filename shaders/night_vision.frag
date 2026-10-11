uniform float xscale=0.0, yscale=0.0, aspect_ratio=1.0;
uniform sampler2D frame_buffer_tex;

in vec2 tc;

void main() {
	vec3 color      = texture(frame_buffer_tex, vec2(tc.x*xscale, tc.y*yscale)).rgb;
	float luminance = (color.r + color.g + color.b)/3.0;
	float gscale    = mix(3.0*luminance, luminance, luminance);
	// add two black bordered circular openings to the left and right
	float ar = 1.7;
	float xl = 0.5/aspect_ratio;
	vec2 xy_scale = vec2(aspect_ratio, 1.0);
	gscale *= clamp(20.0*(0.56 - min(length((tc - vec2(xl, 0.5))*xy_scale), length((tc - vec2((1.0 - xl), 0.5))*xy_scale))), 0.0, 1.0);
	fg_FragColor    = vec4(0.0, gscale, 0.0, 1.0); // all green
}
