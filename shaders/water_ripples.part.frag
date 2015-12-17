uniform sampler2D ripple_tex;

vec3 compute_ripple(in vec2 uv, in float time, in float weight) {
	vec4 ripple        = texture(ripple_tex, uv);
	//vec4 ripple        = textureLod(ripple_tex, uv, 0);
    ripple.yz          = ripple.yz * 2 - 1; // decompress perturbation
    float drop_fract   = fract(ripple.w + time); // apply time shift
    float time_fract   = drop_fract - 1.0 + ripple.x;
    float drop_factor  = clamp((0.2 + weight * 0.8 - drop_fract), 0.0, 1.0);
    float final_factor = drop_factor * ripple.x * sin(clamp(time_fract * 9.0, 0.0, 3.0) * 3.14159);
    return vec3(ripple.yz * final_factor * 0.35, 1.0);
}

vec3 get_ripple_normal(in vec2 uv, in float time, in float rain_intensity) {
	vec4 time_mul = vec4(1.0, 0.85, 0.93, 1.13);
	vec4 time_add = vec4(0.0, 0.2,  0.45, 0.7 );
	vec4 times    = (0.05*time * time_mul + time_add) * 1.6;
	times         = fract(times);
	vec4 weights  = rain_intensity - vec4(0, 0.25, 0.5, 0.75);
	weights       = clamp(4.0*weights, 0.0, 1.0);
	// generate four shifted layer of animated circle
	vec3 ripple1 = compute_ripple(uv + vec2( 0.25, 0.0), times.x, weights.x);
	vec3 ripple2 = compute_ripple(uv + vec2(-0.55, 0.3), times.y, weights.y);
	vec3 ripple3 = compute_ripple(uv + vec2(0.6,  0.85), times.z, weights.z);
	vec3 ripple4 = compute_ripple(uv + vec2(0.5, -0.75), times.w, weights.w);
	// enable one layer by quarter intensity and progressively blend in the current layer
	// compose normal of the four layer based on weights
	vec4 z = mix(vec4(1.0), vec4(ripple1.z, ripple2.z, ripple3.z, ripple4.z), weights);
	return vec3(weights.x*ripple1.xy + weights.y*ripple2.xy + weights.z*ripple3.xy + weights.w*ripple4.xy, z.x*z.y*z.z*z.w);
}

