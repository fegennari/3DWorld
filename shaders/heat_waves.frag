uniform sampler2D frame_buffer_tex;
uniform float time = 0.0;

in vec2 tc;

float hash(vec2 p) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
float hash(vec3 p) {
	float h = dot(p,vec3(127.1,311.7,0.3));	
    return fract(sin(h)*43758.5453123);
}

void main() {
	vec2 pos     = tc;
	//pos.x       += 0.01*sin(20.0*pos.x + 10.0*pos.y + 0.01*time); // drunk mode
	pos.x       += 0.01*sin(124.0*pos.x + 437.0*pos.y - 0.1*time) * (0.5 - abs(tc.x - 0.5)) * (0.5 - abs(tc.y - 0.5));
	vec3 color   = texture(frame_buffer_tex, pos).rgb;
	//color.rgb = color.brg;
	//color *= 10.0;
	//float mag = (color.r + color.g + color.b)/3.0; color = vec3(mag); // convert to grayscale
	//color = mix(color, vec3(hash(vec3(pos, time))), 0.5); // static
	fg_FragColor = vec4(color, 1.0);
}
