in vec4 epos;
in vec3 eye_norm;

void main() {
	fg_FragColor = gl_Color;
	fg_FragColor.a *= pow(min(1.0, 1.2*abs(dot(normalize(eye_norm), normalize(epos.xyz)))), 2); // fade at cone edges (has banding)
	//float d_delta   = log_to_linear_depth(get_depth_at_fragment()) - log_to_linear_depth(gl_FragCoord.z); // in depth_utils.part.frag
	//fg_FragColor.a *= smoothstep(0.0, 1.0, clamp(d_delta/0.1, 0.0, 1.0));
}
