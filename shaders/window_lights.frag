uniform sampler2D tex0;
uniform float min_alpha = 0.0;

in vec2 tc;

float rand(vec2 co) {
   return fract(sin(dot(co.xy,vec2(12.9898,78.233))) * 43758.5453);
}
highp float rand_v2(vec2 co) { // from http://byteblacksmith.com/improvements-to-the-canonical-one-liner-glsl-rand-for-opengl-es-2-0/
    highp float a = 12.9898;
    highp float b = 78.233;
    highp float c = 43758.5453;
    highp float dt= dot(co.xy ,vec2(a,b));
    highp float sn= mod(dt,3.14);
    return fract(sin(sn) * c);
}

void main() {
	vec4 color = texture(tex0, tc);
	if (color.a <= min_alpha) discard;
	// use integer component of window texture as hash to select random lit vs. dark windows
	float rand_val = rand_v2(vec2(floor(tc.s), floor(tc.y)));
	if (rand_val < 0.67) discard; // 67% of windows are dark
	fg_FragColor = color;
}
