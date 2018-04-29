uniform sampler2D tex0;
uniform float min_alpha = 0.0;

in vec2 tc;
in vec4 epos;

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
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard; // transparent part
	// use integer component of window texture as hash to select random lit vs. dark windows
	float rand_val = rand_v2(vec2(floor(tc.s), floor(tc.y)));
	float lit_thresh = fract(5.67*sin(12.3*gl_Color.b)); // use per-building color as hash for lit percentage
	if (rand_val < (0.5 + 0.4*lit_thresh)) discard; // 50% - 90% of windows are dark
	float alpha     = 1.0 - 2.0*(texel.r - 0.5); // mask off window border (white), keep window pane (gray)
	alpha          *= min(1.0, (0.5*length(epos.xyz) - 0.25)); // increase alpha when very close to the camera, as nearby windows don't look as good
	if (alpha < 0.01) discard;
	float rgb_scale = 1.0/max(gl_Color.r, max(gl_Color.g, gl_Color.b)); // normalize largest color component to 1.0
	fg_FragColor    = apply_fog_epos(vec4(rgb_scale*gl_Color.rgb, 1.0), epos);
	fg_FragColor.a *= alpha*texel.a; // keep orig alpha value
}
