uniform float fog_scale   = 0.0;
uniform float fog_density = 1.0;
uniform vec4 fog_color    = vec4(1.0);

// exponential fog
vec4 apply_fog(in vec4 color) {
	float fog = exp(-fog_density*gl_FogFragCoord);
	return vec4(mix(fog_color.rgb, color.rgb, mix(1.0, fog, fog_scale)), color.a);
}
