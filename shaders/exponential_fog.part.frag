uniform float fog_scale = 0.0;

// exponential fog
vec4 apply_fog(in vec4 color) {
	float fog = exp(-gl_Fog.density*gl_FogFragCoord);
	return vec4(mix(gl_Fog.color.rgb, color.rgb, mix(1.0, fog, fog_scale)), color.a);
}
