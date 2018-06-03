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
	float wind_id_y  = floor(tc.t); // window floor/row index
	float max_gsz    = 1.0 + 10.0*fract(6.78*sin(45.6*gl_Color.g)); // max window group size, per-building, 1-10
	max_gsz         *= 0.5 + 0.5*fract(8.91*sin(13.7*wind_id_y)); // modulate 50% by floor
	float group_sz   = max(1.0, ceil(max_gsz*rand_v2(vec2(wind_id_y, 0.001*gl_Color.g)))); // create lit/dark rows of windows
	float wind_id_x  = floor(tc.s/group_sz); // window column index
	float rand_val   = rand_v2(vec2(wind_id_x, wind_id_y));
	float lit_thresh = fract(5.67*sin(12.3*gl_Color.b)); // use per-building color as hash for lit percentage
	lit_thresh      *= 0.5 + 1.0*fract(8.91*sin(13.7*wind_id_y)); // modulate by floor
	float dist       = length(epos.xyz);
	float alpha      = 1.0 - 2.0*(texel.r - 0.5); // mask off window border (white), keep window pane (gray)
	if (dist < 4.0) {alpha *= clamp(3.0*(fract(tc.t) - 0.1), 0.0, 1.0);} // top of window is brighter than bottom
	if (rand_val < (0.5 + 0.4*lit_thresh)) {alpha *= 0.2*rand_val;} // 50% - 90% of windows are darker
	float rand_amt = 0.2*(1.0 - alpha) + 0.15;
	alpha *= min(1.0, (0.5*dist - 0.25)); // increase alpha when very close to the camera, as nearby windows don't look as good
	if (alpha < 0.01) discard;
	float rgb_scale = 1.0/max(gl_Color.r, max(gl_Color.g, gl_Color.b)); // normalize largest color component to 1.0
	fg_FragColor    = apply_fog_epos(vec4(rgb_scale*gl_Color.rgb, 1.0), epos);
	fg_FragColor.a *= alpha*texel.a; // keep orig alpha value
	// add some random per-window color variation
	fg_FragColor.bg *= mix(1.0, fract(13.5*rand_val), rand_amt);
	fg_FragColor.b  *= mix(1.0, fract(17.3*rand_val), rand_amt);
}
